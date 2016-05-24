tc = function(X){t(Conj(X))}

# @export
speclagreg.K.experimental = function(X,Y,freq,method = 0,SXX = NULL){  
  if (method == 2){
    K = freq
    n = dim(X)[1]
    th = 1/sqrt(n)
    for (i in 1:length(K)){
      E = eigen(SXX$operators[i,,])
      K[i] = max(1,sum(abs(E$values) > th))
    }
    
    return(K)    
  }
  
  
  T = length(freq)
  p = 10
  n = dim(X)[1]
  L = floor(n/p)
  Xar = array(0,c(p,L,dim(X)[2]))
  Yar = array(0,c(p,L,dim(Y)[2]))
  FX = array(0,c(L,T,dim(X)[2]))
  FY = array(0,c(L,T,dim(Y)[2]))
  CX = array(0,c(T,dim(X)[2],dim(X)[2]))
  CYX = array(0,c(T,dim(Y)[2],dim(X)[2]))

  # split data into blocks
  for (i in 1:p){
    Xar[i,,] = X[(0:(L-1))*p + i,]
    Yar[i,,] = Y[(0:(L-1))*p + i,]
  }

  # compute FT of each block
  for (i in 1:L){
    FX[i,,] = fourier.transform(Xar[,i,],freq=freq)$operators[,,1] / sqrt(p)
    FY[i,,] = fourier.transform(Yar[,i,],freq=freq)$operators[,,1] / sqrt(p)
  }
  
  # prepare covariances
  # - in case of crossvalidation take 80% of the set for training
  # - in case of our method take the whole dataset for training
  for (i in 1:T){
    if (method == 0){
      tr = 1:(floor(L*0.8))
      ts = (length(tr)+1):(L)
      
      CX[i,,] = t(FX[tr,i,]) %*% Conj(FX[tr,i,])
      CYX[i,,] = t(FY[tr,i,]) %*% Conj(FX[tr,i,])

    }
    else {
      CX[i,,] = t(FX[,i,]) %*% Conj(FX[,i,])
      CYX[i,,] = t(FY[,i,]) %*% Conj(FX[,i,])
    }
  }
  CXfd = freqdom(CX,freq=freq)
  CYXfd = freqdom(CYX,freq=freq)
  #E = freqdom.eigen(CXfd)
  
  K = rep(0,T)
  d = min(dim(X)[2],L-1)
  
  RES = array(0,c(T,d+1))
  for (i in 1:T){
    minMSE = 100000
    for (k in 0:d){
      D = CYX[i,,] %*% pseudoinverse(CX[i,,],K=k)
      FY.est = t(D %*% t(FX[,i,]))

      if (method == 0){
        RES[i,k+1] = MSE(FY[ts,i,],FY.est[ts,])
      }

      if (method == 1){

        res = t(FY.est - FY[,i,])
        cm = rowMeans(res)
        res = res - cm
        cur = abs(mean(diag(tc(res) %*% (res))))/(L - k)^(2) #+ k/L1
        RES[i,k+1] = cur
      }
    }
  }
  K = apply(RES,1,which.min) - 1
  K

}
