unisegclust <- function(Y,CGHo,Kmax){

  P            = CGHo["nblevels"]
  loglik       = c(-Inf,Kmax)
  n.com        = length(Y)
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  mx           = mean(x)
  tol          = 1e-2
  mBIC         = matrix(-Inf,Kmax)
  SSall        = sum((x-mean(x))^2)
  n            = length(x)
  gamma.coef   = lgamma(0.5*(n-P-1)) - lgamma(0.5*(n-1))
 

######   If model selection is performed run P:Kmax and select Kmax
  
  if (CGHo["select"]!="none"){
    select.tmp   = CGHo["select"]
    select(CGHo) = "none"
    for (K in P:Kmax){
      mu        = unisegmean(x,CGHo,K)$mu
      phi       = EMinit(x,as.matrix(mu[,-3]),P,vh=TRUE)
      out.EM    = EMalgo(x,phi,as.matrix(mu[,-3]),P, vh = TRUE)
      mu$mean   = out.EM$phi[apply(out.EM$tau,1,which.max)]
      iter      = 0
      pred      = rep(mu$mean,c(mu$end-mu$begin+1))
      eps       = Inf      
      while ( (eps > tol) & (iter < CGHo@itermax) & (out.EM$empty==0) ){          
        iter     = iter+1
        pred.tmp = pred
        t.est    = ClassiSeg(x,out.EM$phi[1:P],K)$t.est[K,]
        mu[,1:2] = cbind(c(1,t.est[1:(K-1)]+1),t.est)
        out.EM   = EMalgo(x,out.EM$phi,as.matrix(mu[,-3]), P, vh = TRUE)
        mu$mean  = out.EM$phi[apply(out.EM$tau,1,which.max)]
        pred     = rep(mu$mean,c(mu$end-mu$begin+1))
        eps      = max(abs((pred-pred.tmp)/(pred)))
      }
      if (out.EM$empty!=0){
        cat("[unisegclust ERROR]: convergence to a solution with empty levels.","\n");
        stop("[unisegclust ERROR]: try a lower nblevels(CGHo)","\n");
      }
      mu$levels  = apply(out.EM$tau,1,which.max)
      mBIC[K]    = getmBIC(K,out.EM$lvinc,list(aux=mu),CGHo)
    }

    K         = which.max(mBIC)
    mu        = unisegmean(x,CGHo,K)$mu
    phi       = EMinit(x,as.matrix(mu[,-3]),P,vh=TRUE)
    out.EM    = EMalgo(x,phi,as.matrix(mu[,-3]),P, vh = TRUE)
    mu$mean   = out.EM$phi[apply(out.EM$tau,1,which.max)]
    iter      = 0
    pred      = rep(mu$mean,c(mu$end-mu$begin+1))
    eps       = Inf      
    while ( (eps > tol) & (iter < CGHo@itermax) & (out.EM$empty==0) ){          
      iter     = iter+1
      pred.tmp = pred
      t.est    = ClassiSeg(x,out.EM$phi[1:P],K)$t.est[K,]
      mu[,1:2] = cbind(c(1,t.est[1:(K-1)]+1),t.est)
      out.EM   = EMalgo(x,out.EM$phi,as.matrix(mu[,-3]), P, vh = TRUE)
      mu$mean  = out.EM$phi[apply(out.EM$tau,1,which.max)]
      pred     = rep(mu$mean,c(mu$end-mu$begin+1))
      eps      = max(abs((pred-pred.tmp)/(pred)))
    }
    if (out.EM$empty!=0){
      cat("[unisegclust ERROR]: convergence to a solution with empty levels.","\n");
      stop("[unisegclust ERROR]: try a lower nblevels(CGHo)","\n");
    }
    mu$levels    = apply(out.EM$tau,1,which.max)
    select(CGHo) = select.tmp
  } else {
    
######   If no model selection is performed run Kmax
    K         = Kmax
    mu        = unisegmean(x,CGHo,K)$mu
    phi       = EMinit(x,as.matrix(mu[,-3]),P,vh=TRUE)
    out.EM    = EMalgo(x,phi,as.matrix(mu[,-3]),P, vh = TRUE)
    mu$mean   = out.EM$phi[apply(out.EM$tau,1,which.max)]
    iter      = 0
    pred      = rep(mu$mean,c(mu$end-mu$begin+1))
    eps       = Inf      
    while ( (eps > tol) & (iter < CGHo@itermax) & (out.EM$empty==0) ){          
      iter     = iter+1
      pred.tmp = pred
      t.est    = ClassiSeg(x,out.EM$phi[1:P],K)$t.est[K,]
      mu[,1:2] = cbind(c(1,t.est[1:(K-1)]+1),t.est)
      out.EM   = EMalgo(x,out.EM$phi,as.matrix(mu[,-3]), P, vh = TRUE)
      mu$mean  = out.EM$phi[apply(out.EM$tau,1,which.max)]
      pred     = rep(mu$mean,c(mu$end-mu$begin+1))
      eps      = max(abs((pred-pred.tmp)/(pred)))
    }
    if (out.EM$empty!=0){
      cat("[unisegclust ERROR]: convergence to a solution with empty levels.","\n");
      stop("[unisegclust ERROR]: try a lower nblevels(CGHo)","\n");
    }
    mu$levels  = apply(out.EM$tau,1,which.max)
  }
  
######   output
  
  th      = bpwmissing.calls(t.est,present.data,n.com)  
  rupt    = matrix(ncol=2,c(c(1,th[1:length(th)-1]+1),th))    
  mu      = data.frame(begin = rupt[,1],
    end   = rupt[,2],
    mean  =  out.EM$phi[apply(out.EM$tau,1,which.max)],
    levels = mu$levels)
  invisible(list(mu=mu,loglik=out.EM$lvinc))        
}

bpwmissing.calls <- function(t.est,present.data,n.com){
  K = length(t.est)
  for (h in 1:K){
    if (length(which(t.est[h]==0))!=0){
      t.est[h][-which(t.est[h]==0)] = present.data[t.est[h]]
    } 
    else {t.est[h] = present.data[t.est[h]]}
  }
  t.est[K] = n.com
  invisible(t.est)
}



