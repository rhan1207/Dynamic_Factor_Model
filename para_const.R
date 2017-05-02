para_const <- function(X,P,lag){
  
  Z_0 = P$Z_0
  V_0 = P$V_0
  A = P$A
  C = P$C
  Q = P$Q
  R = P$R
  Mx = P$Mx
  Wx = P$Wx
  
  
  ###########################
  # Preparation of the data #
  ###########################
  
  T = dim(X)[1]
  N = dim(X)[2]
  
  # Standaridise x
  
  xNaN = (X - matrix(rep(Mx, T), T, dim(Mx)[2], byrow = TRUE)) / matrix(rep(Wx, T), T, dim(Wx)[2], byrow = TRUE) 
  #source("runKF.R")
  y = t(xNaN)
  Sf = SKF(y,C,R,A,Q,Z_0,V_0)
  Ss = FIS(y,C,R,A,Q,Sf)
  
  Ps = Ss$PmT[,,-1]
  Pf = Sf$PmU[,,-1]
  
  Plag = list()
  Plag[[1]] = Ps
  if (lag > 0){
    for (jk in 1:lag){
      Plag[[jk+1]] = array(0, dim(Ps))
      for (jt in dim(Ps)[3]:(lag+1)){
        As = Pf[,,(jt-jk)] %*% t(A) %*% ginv(A %*% Pf[,,(jt-jk)] %*% t(A) + Q)
        Plag[[jk+1]][,,jt] = As %*% Plag[[jk]][,,jt]
      }
    }
  }
  Zsmooth = Ss$AmT
  Vsmooth = Ss$PmT
  
  Zsmooth = t(Zsmooth)
  x_sm = Zsmooth[2:dim(Zsmooth)[1],] %*% t(C)
  X_sm = matrix(rep(Wx, T), T, dim(Wx)[2], byrow = TRUE) * x_sm + matrix(rep(Mx, T), T, dim(Mx)[2], byrow = TRUE)
  
  Res = list(Plag, Vsmooth, X_sm, Zsmooth[2:dim(Zsmooth)[1],])
  names(Res) = c("Plag", "P", "X_sm", "F")
  return(Res)
  
}