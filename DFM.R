DFM <- function(X, Par){
  #source("runKF.R")
  t = dim(X)[1]  
  N = dim(X)[2]
  
  r = Par$r                      # number of components we want to extract for each block
  p = Par$p
  nQ = Par$nQ
  blocks = Par$blocks
  i_idio = as.logical(c(rep(1, N-nQ), rep(0, nQ)))
  
  #R*Lambda = q; Contraints on the loadings of the quartrly variables
  R_mat = matrix(c(2,-1,0,0,0, 3,0,-1,0,0, 2,0,0,-1,0, 1,0,0,0,-1), nrow=4, ncol=5, byrow = TRUE)
  q = rep(0,4)
  thresh = 1e-5
  max_iter = 5000
  ### Preparation of the data, standarize
  Mx = t(matrix(apply(X, 2, mean, na.rm = TRUE)))
  Wx = t(matrix(apply(X, 2, sd, na.rm = TRUE)))
  xNaN = (X - matrix(rep(Mx, each = t), nrow = t)) / matrix(rep(Wx, each = t), nrow = t)
  
  ### Initial Conditions
  optNaN = "options"
  optNaN$method = 2
  optNaN$k = 3          # parameter for Moving Average
  
  theta = InitCond(xNaN, r, p, blocks, optNaN, R_mat, q, nQ, i_idio) 
  
  A   = theta$A
  C   = theta$C
  Q   = theta$Q
  R   = theta$R
  Z_0 = theta$initZ
  V_0 = theta$initV

  # some auxiliary variables for the iterations
  previous_loglik = -Inf
  num_iter = 0
  LL = -Inf
  converged = 0
  
  # y for the estimation is WITH missing data
  y = t(xNaN)
  
  
  ############################################################################
  # THE EM LOOP
  ############################################################################
  
  #The model can be written as 
  # y = C*Z + e
  # Z = A*Z(t-1) + v
  #where y is N x T, Z is (pr) x T, etc
  
  #remove the leading and ending nans for the estimation
  optNaN$method = 3
  y_est = t(remNaNs(xNaN, optNaN)$X)
  
  decrease = 1:max_iter

  while (num_iter < max_iter && !converged) {
    em_res = EMstep(y_est, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ, i_idio, blocks)
    C = em_res$C_new
    R = em_res$R_new
    A = em_res$A_new
    Q = em_res$Q_new
    Z_0 = em_res$Z_0
    V_0 = em_res$V_0
    loglik = em_res$loglik
    # Checking convergence
    if (num_iter > 2){
      conv_res = em_converged(loglik, previous_loglik, thresh)
      converged = conv_res$converged
      decrease[num_iter+1] = conv_res$decrease      
    }
    
    LL = cbind(LL, loglik)
    previous_loglik = loglik
    num_iter = num_iter + 1
    print(num_iter)
  }
  
  #final run of the Kalman filter
  Zsmooth = runKF(y, A, C, Q, R, Z_0, V_0)
  sm = t(Zsmooth$xsmooth)
  x_sm = sm[2:dim(sm)[1],] %*% t(C)
  X_sm = matrix(rep(Wx, t), t, dim(Wx)[2], byrow = TRUE) * x_sm + matrix(rep(Mx, t), t, dim(Mx)[2], byrow = TRUE) 


  Res = list(X_sm, sm[2:dim(sm)[1], ], C, R, A, Q, Mx, Wx, Z_0, V_0, r, p)
  names(Res) = c("X_sm", "F", "C", "R", "A", "Q", "Mx", "Wx", "Z_0", "V_0", "r", "p")
  
  return(Res)
}


EMstep <- function (y, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ, i_idio, blocks){
  
n = dim(y)[1]
TT = dim(y)[2]

nM = n-nQ
  
pC = dim(R_mat)[2]
ppC = max(p,pC)

n_b = dim(blocks)[2]
# Compute the (expected) sufficient statistics for a single Kalman filter sequence.
# Running the Kalman filter with the current estimates of the parameters
kf_res = runKF(y, A, C, Q, R, Z_0, V_0)

Zsmooth  = kf_res$xsmooth
Vsmooth  = kf_res$Vsmooth
VVsmooth = kf_res$VVsmooth
loglik   = kf_res$loglik

A_new = A
Q_new = Q
V_0_new = V_0

for (i in 1:n_b){
  
  r_i = r[i]
  rp = r_i * p
  rp1 = 0
  if (i > 1){
    rp1 = sum(r[1:(i-1)]) * ppC
  }
  
  A_i = A[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
  Q_i = Q[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
  
  EZZ = Zsmooth[(rp1+1):(rp1+rp), -1, drop = FALSE] %*%  t(Zsmooth[(rp1+1):(rp1+rp), -1, drop = FALSE]) + apply(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),-1, drop = FALSE], c(1,2), sum) 
  EZZ_BB = Zsmooth[(rp1+1):(rp1+rp), -dim(Zsmooth)[2], drop = FALSE] %*%  t(Zsmooth[(rp1+1):(rp1+rp), -dim(Zsmooth)[2],drop = FALSE]) + apply(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),-dim(Zsmooth)[2],drop = FALSE], c(1,2),sum) 
  EZZ_FB = Zsmooth[(rp1+1):(rp1+rp), -1, drop = FALSE] %*%  t(Zsmooth[(rp1+1):(rp1+rp), -dim(Zsmooth)[2],drop = FALSE]) + apply(VVsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),,drop = FALSE],c(1,2),sum) 

  A_i[1:r_i, 1:rp] = EZZ_FB[1:r_i, 1:rp, drop = FALSE] %*% solve(EZZ_BB[1:rp, 1:rp])
  Q_i[1:r_i, 1:r_i] = (EZZ[1:r_i, 1:r_i] - A_i[1:r_i,1:rp,drop = FALSE] %*% t(EZZ_FB[1:r_i,1:rp, drop = FALSE])) / TT
  
  A_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] = A_i
  Q_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] = Q_i
  V_0_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] = Vsmooth[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC),1]
}

rp1 = sum(r) * ppC
niM = sum(i_idio[1:nM])

#idiosyncratic
EZZ = diag(diag(Zsmooth[(rp1+1):dim(Zsmooth)[1], -1, drop = FALSE] %*%  t(Zsmooth[(rp1+1):dim(Zsmooth)[1], -1, drop = FALSE]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],-1, drop = FALSE], c(1,2), sum))) 
EZZ_BB = diag(diag(Zsmooth[(rp1+1):dim(Zsmooth)[1], -dim(Zsmooth)[2], drop = FALSE] %*%  t(Zsmooth[(rp1+1):dim(Zsmooth)[1], -dim(Zsmooth)[2],drop = FALSE]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],-dim(Vsmooth)[3],drop = FALSE], c(1,2),sum))) 
EZZ_FB = diag(diag(Zsmooth[(rp1+1):dim(Zsmooth)[1], -1, drop = FALSE] %*%  t(Zsmooth[(rp1+1):dim(Zsmooth)[1], -dim(Zsmooth)[2],drop = FALSE]))) + diag(diag(apply(VVsmooth[(rp1+1):dim(VVsmooth)[1],(rp1+1):dim(VVsmooth)[2],,drop = FALSE],c(1,2),sum))) 


A_i = EZZ_FB %*% diag(1/diag(EZZ_BB))
Q_i = (EZZ - A_i %*% t(EZZ_FB)) / TT

A_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] = A_i[1:niM, 1:niM]
Q_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] = Q_i[1:niM, 1:niM]
V_0_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] = diag(diag(Vsmooth[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM),1]))

Z_0 = Zsmooth[,1]

nanY = is.na(y)
y[nanY] = 0


#LOADINGS 
C_new = C
#Blocks
bl = unique(blocks)
temp_fix = c(3,1,4,2)
bl = bl[temp_fix, ]    # have to find a better way to handle this...
n_bl = dim(bl)[1]

for (i in 1:n_b){
  if (i == 1) {
    bl_idxQ = matrix(rep(bl[,i],r[i]*ppC), ncol = r[i]*ppC, byrow = FALSE)
    bl_idxM = cbind(matrix(rep(bl[,i],r[i]), ncol = r[i], byrow = FALSE), matrix(0, n_bl, r[i] * (ppC-1)))
    R_con = kronecker(R_mat, diag(r[i]))
    q_con = matrix(0, dim(R_mat)[1]*r[i], 1)
  }
  else{
    bl_idxQ = cbind(bl_idxQ, matrix(rep(bl[,i],r[i]*ppC), ncol = r[i]*ppC, byrow = FALSE))
    bl_idxM = cbind(bl_idxM, cbind(matrix(rep(bl[,i],r[i]), ncol = r[i], byrow = FALSE), matrix(0, n_bl, r[i] * (ppC-1))))  
    #rr = dim(kronecker(R_mat, diag(r[i])))[1]
    #cc = dim(kronecker(R_mat, diag(r[i])))[2]
    R_con = adiag(R_con, kronecker(R_mat, diag(r[i])))
    #R_con = rbind(cbind(R_con, matrix(0, rr*(i-1), cc)), cbind(matrix(0, rr, cc*(i-1)), kronecker(R_mat, diag(r[i]))))
    q_con = rbind(q_con, matrix(0, dim(R_mat)[1]*r[i], 1))
  }    
}  
  

#bl_idxM = as.logical(bl_idxM)
#bl_idxQ = as.logical(bl_idxQ)

#idio
i_idio_M = i_idio[1:nM]
n_idio_M = length(which(i_idio_M))
c_i_idio = cumsum(i_idio)


for (i in 1:n_bl){
  
  bl_i = bl[i,]
  rs = sum(r[which(bl_i == 1)])
  idx_i = which(apply(blocks, 1, identical, bl_i))
  
  
  #MONTHLY
  idx_iM = idx_i[which(idx_i <  nM + 1)]
  n_i = length(idx_iM)
  
  denom = matrix(0, n_i*rs, n_i*rs)
  nom = matrix(0, n_i, rs)
  
  i_idio_i = i_idio_M[idx_iM]
  i_idio_ii = c_i_idio[idx_iM]
  i_idio_ii = i_idio_ii[i_idio_i]
  temp_idx = which(bl_idxM[i,] == 1)
   for (t in 1:TT){
     nanYt = diag(!nanY[idx_iM,t])
     denom = denom + kronecker(Zsmooth[temp_idx, t+1, drop = FALSE] %*% t(Zsmooth[temp_idx, t+1, drop = FALSE]) +
             Vsmooth[temp_idx, temp_idx, t+1], nanYt)
     
     ###### ?????? shouldn't be + ###
     nom = nom + y[idx_iM,t,drop = FALSE] %*% t(Zsmooth[temp_idx,t+1,drop = FALSE]) - 
           nanYt[,i_idio_i, drop = FALSE] %*% (Zsmooth[rp1+i_idio_ii, t+1,drop = FALSE] %*% t(Zsmooth[temp_idx, t+1, drop = FALSE])
           + Vsmooth[rp1+i_idio_ii, temp_idx,t+1])
   }
   vec_C = solve(denom)%*%vecFun(nom)
   C_new[idx_iM, which(bl_idxM[i,] == 1)] = matrix(vec_C, n_i, rs)
  
  # QUARTERLY
  idx_iQ = idx_i[which(idx_i > nM)]
  rps = rs * ppC

  R_con_i = R_con[, which(bl_idxQ[i,] == 1)]
  q_con_i = q_con
  
  no_c = apply(R_con_i, 1, any)
  R_con_i = R_con_i[which(no_c),]
  q_con_i = q_con_i[which(no_c),]

  for (j in idx_iQ){
    denom = matrix(0, rps, rps)
    nom = matrix(0, 1, rps)
    idx_jQ = j - nM
    i_idio_jQ = (rp1 + n_idio_M + 5*(idx_jQ-1) + 1) : (rp1+n_idio_M + 5*idx_jQ)
    V_0_new[i_idio_jQ, i_idio_jQ] = Vsmooth[i_idio_jQ, i_idio_jQ, 1]
    A_new[i_idio_jQ[1], i_idio_jQ[1]] = A_i[i_idio_jQ[1]-rp1, i_idio_jQ[1] - rp1]
    Q_new[i_idio_jQ[1], i_idio_jQ[1]] = Q_i[i_idio_jQ[1]-rp1, i_idio_jQ[1] - rp1]
     for (t in 1:TT) {
       nanYt = as.integer(!nanY[j,t])
       denom = denom + kronecker(Zsmooth[which(bl_idxQ[i,] == 1), t+1, drop = FALSE] %*% 
                                   t(Zsmooth[which(bl_idxQ[i,] == 1),t+1, drop = FALSE]) +
                                   Vsmooth[which(bl_idxQ[i,] == 1), which(bl_idxQ[i,] == 1), t+1], nanYt)
      nom = nom + y[j,t,drop = FALSE] %*% t(Zsmooth[which(bl_idxQ[i,] == 1), t+1, drop = FALSE])
      nom = nom - nanYt %*% (matrix(c(1,2,3,2,1), 1, 5) %*% Zsmooth[i_idio_jQ, t+1, drop = FALSE] %*%
            t(Zsmooth[which(bl_idxQ[i,] == 1), t+1, drop = FALSE]) + 
            matrix(c(1,2,3,2,1), 1, 5) %*% Vsmooth[i_idio_jQ, which(bl_idxQ[i,] == 1), t+1])
      
    }
    C_i = solve(denom) %*% t(nom)
    C_i_constr = C_i - solve(denom) %*% t(R_con_i) %*% solve(R_con_i %*% solve(denom) %*% t(R_con_i)) %*%
                (R_con_i %*% C_i - q_con_i)
    C_new[j, which(bl_idxQ[i,] == 1)] = C_i_constr
  }
  
}

R_new = matrix(0, n, n)
for (t in 1:TT){
  nanYt = diag(!nanY[,t])
  R_new = R_new + (y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) %*% t(y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) +
    nanYt %*% C_new %*% Vsmooth[,,t+1] %*% t(C_new) %*% nanYt + 
    (diag(n)-nanYt) %*% R %*% (diag(n) - nanYt)    
}

R_new = R_new / TT
RR = diag(R_new)
RR[i_idio_M] = 1e-04
if (nQ > 1){RR[(nM+1):length(RR)] = 1e-04}
R_new = diag(RR)

res = list(C_new, R_new, A_new, Q_new, Z_0, V_0, loglik)
names(res) = c("C_new", "R_new", "A_new", "Q_new", "Z_0", "V_0", "loglik")
return(res)
}



em_converged <- function(loglik, previous_loglik, threshold, check_increased){
# EM_CONVERGED has EM converged?
# converged if the slope of the log-likelihood under 'threshold'   
# i.e. |f(t) - f(t-1)| / avg < threshold,
# where avg = (|f(t)| + |f(t-1)|) / 2 and f(t) is log lik at iteration t.
# 'threshold' defaults to 1e-4
#
# This stopping criterion is from Numerical Recipes in C p423
# 
# If we are doing MAP estimation (using priors), the likelihood can decrease  
# even though the mode of the posterior is increasing
  
if (missing(threshold)){ threshold = 1e-4}
if (missing(check_increased)){ check_increased = 1}

converged = 0
decrease = 0

if (check_increased){
  if ((loglik - previous_loglik) < -1e-3){
    sprintf("******likelihood decreased from %6.4f to %6.4f!\n", previous_loglik, loglik)
    decrease = 1
  } 
}
  
delta_loglik = abs(loglik - previous_loglik)
avg_loglik = (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
if (delta_loglik / avg_loglik < threshold) {converged  = 1}
  
res = list(converged, decrease)
names(res) = c("converged", "decrease")
return(res)  
}

InitCond <- function(x, r, p, blocks, optNaN, Rcon, q, NQ, i_idio){
pC = dim(Rcon)[2]          # Rcon is the restriction condition
ppC = max(p,pC)
n_b = dim(blocks)[2]       # number of blocks or factors


res <- remNaNs(x,optNaN)
xBal = res$X
indNaN = res$indNaN

t = dim(xBal)[1]
N = dim(xBal)[2]
NM = N - NQ

xNaN = xBal
xNaN[indNaN] = NA                       # reconstruct data with nans

A = matrix(, nrow = 0, ncol = 0)         # State Eq Coeffcieints
Q = matrix(, nrow = 0, ncol = 0)         # State Eq Cov for factors
initV = matrix(, nrow = 0, ncol = 0)     # State Eq Cov for errors


res = xBal                               
resNaN = xNaN
indNaN[1:pC-1, ] = TRUE
for (i in 1:n_b){
  # r[i] determines the number of PC used in the block
  r_i = r[i]     
  
  ######################################################
  ################ Observation equation ################
  ######################################################
  
  C_i = matrix(0, N, r_i * ppC)
  #Selection Matrix W & W^Q
  idx_i = which(blocks[,i] == 1)
  idx_iM = idx_i[which(idx_i<NM+1)]         # valid monthly data
  idx_iQ = idx_i[which(idx_i>NM)]           # valid quarterly data
  decompose = eigen(cov(res[,idx_iM]))      # extract PC for monthly data
  d = decompose$values[1:r_i]
  v = matrix(decompose$vectors[,1:r_i], nrow = length(idx_iM), ncol = r_i)
  
  ###### Use first PC as the inital monthly coefficients #######
  
  C_i[idx_iM, 1:r_i] = v
  
  
  #Multiply original data with the loadings to get the commonfactors
  #EQ(13) in BGR paper
  f = as.matrix(res[,idx_iM]) %*% v  
  # initilize F with kk = 0
  # F contains p+1 or pC columns of lag factors 
  F = f[pC:dim(f)[1],]
  for (kk in 1:(max(p+1, pC)-1)){
     F = cbind(F, f[(pC-kk):(dim(f)[1]-kk),])
  }
  Rcon_i = kronecker(Rcon, diag(r_i))
  q_i = kronecker(q, matrix(0, r_i, 1))
  
  ff = F[,1:(r_i*pC)]   # ff is just F when r_i = 1; five lags of first component 
  
  
  ####### Compute initial quarterly coefficeints #######
  
  for (j in idx_iQ){
    #check how many NAs - if too many, take the sploned data instead...
    xx_j = resNaN[pC:dim(resNaN)[1],j, drop = FALSE]
    if (sum(!is.na(xx_j)) < (dim(ff)[2] + 2)){
      xx_j = res[pC:dim(resNaN)[1],j, drop = FALSE]
    }
    ##monthly factors correspond to quarterly data 
    ff_j = ff[!is.na(xx_j),]
    xx_j = xx_j[!is.na(xx_j)]
    ## map monthly loadings to quarterly data EQ(13) -- D^(-1) matrix 
    iff_j = solve(t(ff_j) %*% ff_j)
    ##Unrestricted quarterly loadings
    Cc = iff_j %*% t(ff_j) %*% xx_j
    ## adding restriction Rcon_i to construct quarterly loadings 
    Cc = Cc - iff_j %*% t(Rcon_i) %*% solve(Rcon_i %*% iff_j %*% t(Rcon_i)) %*% (Rcon_i %*% Cc - q_i)
    C_i[j, 1:(pC*r_i)] = t(Cc) 
  }
  
  
  
  ######## Compute the residual of the Obs Eqs ############
  
  ff = rbind(matrix(0, pC-1, pC*r_i), ff)
  res = res - ff%*%t(C_i)
  resNaN = res
  resNaN[indNaN] = NA
  
  ######## Store Observation Eq Coefficients #######
  if (i == 1){C = C_i}
  else{C = cbind(C, C_i)}
  
  #########################################################
  ################ Transition equation ####################
  #########################################################
  z = F[, 1:r_i]
  Z = F[, (r_i+1):(r_i*(p+1))]
  
  ####### use regression solution as initial coefficients of states #####
  A_temp = solve(t(Z)%*%Z)%*%t(Z)%*%z
  A_i = matrix(0, r_i*ppC, r_i*ppC)

  A_i = rbind(cbind(t(A_temp), matrix(0, r_i, r_i*(ppC-p))), cbind(diag(r_i*(ppC-1)), matrix(0, r_i*(ppC-1), r_i)))

  ######## Covariance Matrix for Residuals ########
  Q_i = matrix(0, ppC*r_i, ppC*r_i)
  e = matrix(z - Z%*%A_temp, ncol = r_i)
  Q_i[1:r_i, 1:r_i] = cov(e)
  
  ######## MSE of States -- factors #########
  initV_i = matrix(solve(diag((r_i*ppC)^2) - kronecker(A_i, A_i))%*%matrix(Q_i, dim(Q_i)[1]*dim(Q_i)[2], 1), r_i*ppC, r_i*ppC)
  A = adiag(A, A_i)
  Q = adiag(Q, Q_i)
  initV = adiag(initV, initV_i)
}


####### Complete Obs and State Eqs here #######
browser()
#### R is covariance of errors in obs eqs, which is calibrated to 1e-04
R <- matrix(0,N,N)
diag(R) = apply(resNaN,2, var,na.rm = TRUE)    ### ??? recalibarated to a very small number
eyeN = diag(N)[,i_idio]
C = cbind(C, eyeN)
C = cbind(C, rbind(matrix(0, NM, ppC*NQ), kronecker(diag(NQ), matrix(c(1,2,3,2,1), 1, 5)))) # it is 5(ppC) in matlab code


ii_idio = which(i_idio)
n_idio = length(ii_idio)
BM = matrix(0, n_idio, n_idio)      
SM = matrix(0, n_idio, n_idio)


##### Compute the coefficients and Coveriance of residuals ######
for (i in 1:n_idio){
  R[ii_idio[i], ii_idio[i]] = 1e-04
  res_i = resNaN[, ii_idio[i]]
  leadZero = t(1:t) == cumsum(is.na(res_i))
  endZero = rev(t(1:t) == cumsum(is.na(rev(res_i))))
  
  exclude_zeros = leadZero|endZero
  res_i = res[,ii_idio[i]]
  res_i = res_i[!exclude_zeros]
 
  len = length(res_i)
  res_i = matrix(res_i, len, 1)
  
  ###### coefs and variance of states -- monthly errors ######
  BM[i,i] = solve(t(res_i[1:(len-1),])%*%res_i[1:(len-1),])%*%t(res_i[1:(len-1),])%*%res_i[2:len,]
  SM[i,i] = cov(matrix(res_i[2:len,] - res_i[1:(len-1),]*BM[i,i], dim(res_i)[1] - 1, 1))
}

######## MSE of States -- monthly errors ##########
initViM = diag(1/diag(diag(dim(BM)[1]) - BM^2))*SM

Rdiag = diag(R)

if (NQ > 0){
  sig_e = Rdiag[(NM+1):N]/19
  Rdiag[(NM+1):N] = 1e-04
}
R = diag(Rdiag)

# alpha for quarterly errors
rho0 = 0.1  
if (NQ == 0) {
  initViQ = matrix(0, 0, 0)
  BQ = matrix(0, 0, 0)
  SQ = matrix(0, 0, 0)
}
else{
  ###### coefs and variance of states -- quarterly errors ######
  BQ = kronecker(diag(NQ), rbind(cbind(rho0, matrix(0,1,4)), cbind(diag(4), matrix(0, 4, 1))))
  temp = matrix(0, 5, 5)
  temp[1,1] = 1
  SQ = kronecker(diag((1-rho0^2)*sig_e),temp)
  
  ######## MSE of States -- quarterly errors ##########
  initViQ = matrix(solve(diag((5*NQ)^2) - kronecker(BQ,BQ))%*%matrix(SQ, dim(SQ)[1]*dim(SQ)[1], 1), 5*NQ, 5*NQ)
}

A = adiag(A, BM, BQ)
Q = adiag(Q, SM, SQ)

initZ = matrix(0, dim(A)[1], 1)
initV = adiag(initV, initViM, initViQ)


res = list(A, C, Q, R, initZ, initV)
names(res) = c("A", "C", "Q", "R", "initZ", "initV")
# A state equation coefficeints: dimension = (Sum_i(r_i*ppC) + NM + 5*NQ) * (Sum_i(r_i*ppC) + NM + 5*NQ)
# Q state equation covariance matrix: dimension, same as A
# C observation equation coefficients: dimension = N * (Sum_i(r_i*ppC) + NM + 5*NQ)
# R variance of series: dimension = N * N
# initV MSE of Kalman

return(res)
}











































