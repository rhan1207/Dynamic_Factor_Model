News_DFM <- function(X_old, X_new, Q, t_fcst, v_news){
  
  r = dim(Q$C)[2]
  T = dim(X_new)[1]
  N = dim(X_new)[2]
  singlenews = matrix(0, 1, N)
  IS_fore = matrix()
  y_old = matrix(0, length(t_fcst), length(v_news))
  y_new = matrix(0, length(t_fcst), length(v_news))
  
  
  if (! is.na(X_new[t_fcst, v_news])){
    
    Res_old = para_const(X_old, Q, 0)
    
    for (i in 1:length(v_news)){
       temp = X_new[t_fcst, v_news[i]] - Res_old$X_sm[t_fcst, v_news[i]]
       singlenews[, v_news[i]] = temp
       y_old[1:length(t_fcst),i] = Res_old$X_sm[t_fcst, v_news[i]]    
       y_new[1:length(t_fcst),i] = X_new[t_fcst, v_news[i]]
    }
    
    actual = matrix()
    fore = matrix()
    filt = matrix()
    t_miss = matrix()
    v_miss = matrix()
    innov = matrix()
    IS_fore = matrix()
  }
  
  else{
  
    Mx = Q$Mx  
    Wx = Q$Wx
  
    miss_old = is.na(X_old)
    miss_new = is.na(X_new)
    temp = miss_old - miss_new
    
    t_miss = which(temp == 1, TRUE)[,1]
    v_miss = which(temp == 1, TRUE)[,2]
    
    if (length(v_miss) == 0){
      Res_old = para_const(X_old, Q, 0)
      Res_new = para_const(X_new, Q, 0)
      
      y_old = Res_old$X_sm[t_fcst, v_news]
      y_new = y_old
      
      actual = matrix()
      fore = matrix()
      filt = matrix()
      t_miss = matrix()
      v_miss = matrix()
      innov = matrix()
      IS_fore = matrix()
      
      res = list(y_old, y_new, singlenews, actual, fore, filt, t_miss, v_miss, innov, IS_fore)
      names(res) = c("y_old", "y_new", "singlenews", "actual", "fore", "filt", "t_miss", "v_miss", "innov", "IS_fore")
      return(res)
    }
    

    lag = t_fcst - t_miss
    
    k = max(abs(lag), max(lag) - min(lag))
    
    C = Q$C
    R = t(Q$R)
    
    n_news = length(lag)         # this could be problematic when t_fcst = multiple days

    Res_old = para_const(X_old, Q, k)
    
    Plag = Res_old$Plag

    Res_new = para_const(X_new, Q, 0)
    y_old = Res_old$X_sm[t_fcst, v_news] 
    y_new = Res_new$X_sm[t_fcst, v_news]
    
    #P = Res_old$P[,,2:dim(P)[3]]
    
    
    for (i in 1:length(lag)){
      h = abs(t_fcst-t_miss[i])
      m = max(t_miss[i], t_fcst)
      
      if(t_miss[i] > t_fcst){
        Pp=Plag[[h+1]][,,m]
      }
      else{
        Pp =t(Plag[[h+1]][,,m])
      }
      p_temp = Pp %*% C[v_miss[i], 1:r]
      
      if (i == 1) {P1 = p_temp}
      else {P1 = cbind(P1, p_temp)}
    }
    
    
    innov = vector()
    for (i in 1:length(t_miss)) {
      
      X_new_norm = (X_new[t_miss[i], v_miss[i]] - Mx[v_miss[i]])/Wx[v_miss[i]]
      X_sm_norm = (Res_old$X_sm[t_miss[i], v_miss[i]] - Mx[v_miss[i]])/Wx[v_miss[i]]  
      innov[i] = X_new_norm - X_sm_norm
      
    }
    
    w_max = max(v_miss)
    WW = matrix(0, w_max, w_max)
    for (i in 1:length(lag)){
        for (j in 1:length(lag)){
          h = abs(lag[i] - lag[j])
          m = max(t_miss[i], t_miss[j])
          
          if (t_miss[j] > t_miss[i]){
            Pp = Plag[[h+1]][,,m]
          }
          else{
            Pp = t(Plag[[h+1]][,,m])
          }
          if ((v_miss[i] == v_miss[j]) & (t_miss[i] != t_miss[j])){
            WW[v_miss[i], v_miss[j]]  = 0
          }
          else{
            WW[v_miss[i], v_miss[j]]  = R[v_miss[i], v_miss[j]]
          }
          p_temp = C[v_miss[i], 1:r, drop = FALSE] %*% Pp %*% t(C[v_miss[j], 1:r, drop = FALSE]) + WW[v_miss[i], v_miss[j]]
          if (j == 1) {p2 = p_temp}
          else {p2 = cbind(p2, p_temp)}
        }
        if (i == 1) {P2 = p2}
        else {P2 = rbind(P2, p2)}
    }
    
    
    rm(temp)
    
    #totnews = array(0, c(1, length(v_news)))
    temp = array(0, c(1, length(innov), length(v_news)))
    gain = array(0, c(1, length(innov), length(v_news)))    # dimension correct???
    
    # loop on v_news
    for (i in 1:length(v_news)){
      #totnews[1,i] 
      temp[1,,i] = Wx[1,v_news[i],drop = FALSE] %*% C[v_news[i], 1:r] %*% P1 %*% solve(P2) * innov
      gain[1,,i] = Wx[1,v_news[i],drop = FALSE] %*% C[v_news[i], 1:r] %*% P1 %*% solve(P2)
    }
    
    
    
    singlenews = array(0, c(max(t_miss-min(t_miss)+1), N, length(v_news)))
    
    actual = matrix(0, N, 1)
    fore = matrix(0, N, 1)
    filt = array(0, c(N, 1, length(v_news)))
    
    for (i in 1:length(innov)){
      
      actual[v_miss[i], 1] = X_new[t_miss[i], v_miss[i]]
      fore[v_miss[i], 1]  = Res_old$X_sm[t_miss[i], v_miss[i]]
      
      for (j in 1:length(v_news)){
        singlenews[(t_miss[i] - min(t_miss) + 1), v_miss[i], j] = temp[1,i,j]
        filt[v_miss[i],,j] = gain[,i,j] / Wx[v_miss[i]]
      }
    }
    singlenews = colSums(singlenews[,,1,drop = FALSE])
    
  }
  res = list(y_old, y_new, singlenews, actual, fore, filt, t_miss, v_miss, innov, IS_fore)
  names(res) = c("y_old", "y_new", "singlenews", "actual", "fore", "filt", "t_miss", "v_miss", "innov", "IS_fore")
  return(res)
}