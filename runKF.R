





runKF <- function(y, A, C, Q, R, x_0, Sig_0){
  S = SKF(y, C, R, A, Q, x_0, Sig_0)
  SS = FIS(y, C, R, A, Q, S)
  res = list(SS$AmT, SS$PmT, SS$PmT_1, SS$loglik)
  names(res) = c("xsmooth", "Vsmooth", "VVsmooth", "loglik")
  return(res)

}




SKF <-function(Y, Z, R, TT, Q, A_0, P_0){
########################################################################  
#
# Kalman filter for stationary systems with time-varying
# system matrices and missing data  
# 
# The model is       y_t      = Z * a_t + eps_t
#                    a_t+1    = T * a_t + u_t
#
########################################################################
#
# INPUT 
#        Y           Data                               (nobs   x n) 
# OUTPUT
#        S.Am        Predicted state vector  A_t|t-1    (nobs   x m)
#        S.AmU       Filtered state vectort  A_t|t      (nobs+1 x m)
#        S.Pm        Predicted covariance of A_t|t-1    (nobs   x m x m)
#        S.PmU       Filtered covariance of  A_t|t      (nobs+1 x m x m)
#        S.loglik    Value of likelihood function 
#
########################################################################

# Output  structure & dimensions

n = dim(Z)[1]
m = dim(Z)[2]
nobs = dim(Y)[2]

Am  = matrix(NA, m, nobs)
AmU = matrix(NA, m, nobs + 1)

Pm  = array(NA, dim=c(m, m, nobs))
PmU = array(NA, dim=c(m, m, nobs + 1))
loglik = 0


################################ SKF ####################################
Au = A_0        # A_0|0 
Pu = P_0        # P_0|0

AmU[, 1]  = Au
PmU[,,1]  = Pu


# The equations used in the following loop can be found in Chapter 13 
# of Time Series Analysis by Hamilton, annotated below
for (t in 1 : nobs){
  # A = A_t|t-1
  # P = P_t|t-1

  A = TT %*% Au                     # [13.2.17]
  P = TT %*% Pu %*% t(TT) + Q       # [13.2.21]
  P = 0.5 * (P + t(P))
  
  # handling missing data here
  newData = MissData(Y[, t], Z, R)  
  y_t = newData$y
  Z_t = newData$C
  R_t = newData$R
  L_t = newData$L
  
  #[13.2.15]
  #[13.2.16]
  if (length(y_t) == 0) {
    Au = A
    Pu = P
  }
  else{
    PZ = P%*%t(Z_t)
    iF = solve(Z_t%*%PZ + R_t)
    PZF = PZ%*%iF
    
    V = y_t - Z_t %*% A
    Au = A + PZF %*% V
    Pu = P - PZF %*%t(PZ)
    Pu = 0.5 * (Pu + t(Pu))
    loglik = loglik + 0.5 * (log(det(iF)) - t(V)%*%iF%*%V)
  }
  
  Am[,t] = A
  Pm[,,t] = P
  
  # Au = A_t|t    
  # Pu = P_t|t
  
  AmU[,t+1] = Au
  PmU[,,t+1] = Pu
}

  if (length(y_t) == 0){
    KZ = matrix(0, m, m)
  }
  else{
    KZ = PZF %*% Z_t
  }
res = list(Am, Pm, AmU, PmU, loglik, KZ)
names(res) = c("Am", "Pm", "AmU", "PmU", "loglik", "KZ")
return(res)
}

##########################################################################
FIS <- function(Y, Z, R, TT, Q, S){
##########################################################################
# Fixed interval smoother (Harvey 1989, p. 154)
# FIS returns the smoothed state vector AmT and its covariance matrix PmT
# Use this in conjunction with function SKF
##########################################################################
# INPUT
#        Y             Data                         (nobs x n)
#        S Estimates from Kalman filter SKF
#        S.Am       :  Estimates     a_t|t-1       (nobs x m)
#        S.Pm       :  P_t|t-1 = Cov(a_t|t-1)      (nobs x m)
#        S.AmU      :  Estimates     a_t|t         (nobs x m x m)
#        S.PmU      :  P_t|t   = Cov(a_t|t)        (nobs x m x m)
# OUTPUT
#        S Smoothed estimates added to above
#        S.AmT      : Estimates      a_t|T         (nobs x m)
#        S.PmT      : P_t|T   =  cov(a_t|T)        (nobs x m x m)
#        S.PmT_1    : Cov(a_ta_t-1|T)
#        where m is the dim of state vector and t = 1...T is time
##########################################################################

m    =  dim(S$Am)[1]
nobs =  dim(S$Am)[2]

AmT = matrix(0, m, nobs+1)
PmT = array(0, dim=c(m, m, nobs+1))
PmT_1 = array(0, dim=c(m, m, nobs))
AmT[,nobs+1] = drop(S$AmU[,nobs+1])
PmT[,,nobs+1] = drop(S$PmU[,,nobs+1])
PmT_1[,,nobs] = (diag(m) - S$KZ) %*% TT %*% drop(S$PmU[,,nobs])
J_2 = drop(S$PmU[,,nobs]) %*% t(TT) %*% ginv(drop(S$Pm[,,nobs]))            #[13.6.11]

for (t in nobs:1){
  
  PmU = drop(S$PmU[,,t])  
  Pm1 = drop(S$Pm[,,t])
  P_T = drop(PmT[,,t+1])
  P_T1= drop(PmT_1[,,t])
  
  J_1 = J_2                                                          
  
  AmT[, t] = S$AmU[,t] + J_1 %*% (AmT[,t+1] - TT %*% S$AmU[,t])             #[13.6.16]
  PmT[,,t] = PmU + J_1 %*% (P_T - Pm1) %*% t(J_1)                           #[13.6.20] 
  
  if (t > 1){
    J_2 = drop(S$PmU[,,t-1]) %*% t(TT) %*% ginv(drop(S$Pm[,,t-1]))
    PmT_1[,,t-1] = PmU %*% t(J_2) + J_1 %*% (P_T1 - TT %*% PmU) %*% t(J_2) 
  }
 }

res = list(S$Am, S$Pm, S$AmU, S$PmU, S$loglik, S$KZ,AmT,PmT,PmT_1)
names(res) = c("Am", "Pm", "AmU", "PmU", "loglik", "KZ","AmT","PmT", "PmT_1")
return(res)

}

MissData <- function(y,C,R){
##########################################################################
# PROC missdata
# PURPOSE: eliminates the rows in y & matrices Z, G that correspond to 
#          missing data (NA) in y
# INPUT:   y        vector of observations at time t               (n x 1)  
#          S        KF system matrices  must contain Z and G  
#
# OUTPUT:  y        vector of observations    (reduced)            (# x 1)
#          Z G      KF system matrices        (reduced)            (# x ?)
#          L        To restore standard dimensions                 (n x #)
#                   where # is the non-nan data in y
#
##########################################################################

ix = !is.na(y)
e = diag(length(y)[1])
L = e[,ix, drop = FALSE]
y = y[ix,drop = FALSE]
C = C[ix, ,drop = FALSE]
R = R[ix,ix,drop = FALSE]

res = list(y, C, R, L)
names(res) = c("y", "C", "R", "L")
return(res)
}







