pweight <- function(Z,X,r,wt,program=3){ # Z:treat, X:X, r:residual (from (g)lm with weight), wt:weight
  
  N = length(r) # number of observations
  
  # combine Z, X and intercept
  ZX1 <- cbind(Z,X,rep(1,len=N))
  
  k = ncol(ZX1)   # number of variables
  
  # create blank matrices
  XWeeWX = XWX = matrix(0,k,k)
  
  # calculate meat of sandwich
  
  if (program==3){
    r.wt = r*wt
    XWeeWX.vec = apply(r.wt*ZX1,1,tcrossprod)
    XWeeWX = matrix(as.numeric(apply(XWeeWX.vec,1,sum)),k,k)
  }
  
  if (program==1){
    XWeeWX.all = array(NA, c(k,k,N))
    r.wt.2 = (r*wt)^2
    for (i in 1:N){
      XWeeWX.all[,,i] = r.wt.2[i]^2*tcrossprod(ZX1[i,])
    }
    XWeeWX = apply(XWeeWX.all,c(1,2),sum)
  }
  
  if (program==2){
    XWeeWX.all = vector("list",N)
    r.wt.2 = (r*wt)^2
    for (i in 1:N){
      XWeeWX.all[[i]] = r.wt.2[i]*tcrossprod(ZX1[i,])
    }
    XWeeWX = Reduce("+",XWeeWX.all)
  }

  if (program==4){
    XWeeWX = matrix(as.numeric(apply(apply(r*wt*ZX1,1,tcrossprod),1,sum)),k,k)
  }
  
  # calculate buns of sandwich
  XWX = t(ZX1)%*%diag(wt,N,N)%*%ZX1
  inv.XWX = solve(XWX)
  
  # make sandwitch (df is n-k assuming homoskedasticity and normal distribution)
  V = N/(N-k)*inv.XWX%*%XWeeWX%*%inv.XWX  
  
  # return SE
  sqrt(V[1,1])
}
