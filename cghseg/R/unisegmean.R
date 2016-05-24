unisegmean <- function(Y,CGHo,Kmax){
  
  n.com        = length(Y)
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  out          = segmeanCO(x,Kmax)
  loglik       = -(length(x)/2)*(log(2*pi*out$J.est/(length(x)))+1)
  
  if (CGHo["select"]=="none"){
    Kselect = Kmax
  } else {    
    mBIC = sapply(1:Kmax,FUN=function(K){
      th      = out$t.est[K,1:K]
      rupt    = matrix(ncol=2,c(c(1,th[1:K-1]+1),th))
	  resmean = meanRuptR_c(Y, rupt[,2], K)
      mu      = data.frame(begin = rupt[,1],
        end   = rupt[,2],
        mean  = resmean)    
      mu      = list(aux=mu)  
      getmBIC(K,loglik[K],mu,CGHo)   
	})
    Kselect = which.max(mBIC)
  }
  
  t.est   = bpwmissing(out$t.est,present.data,n.com)
  th      = t.est[Kselect,1:Kselect]
  rupt    = matrix(ncol=2,c(c(1,th[1:Kselect-1]+1),th)) 
  resmean = meanRuptR_c(Y, rupt[,2], Kselect)
  mu      = data.frame(begin = rupt[,1],
    end   = rupt[,2],
    mean  = resmean)     
  invisible(list(mu=mu,loglik=loglik,t.est=t.est))
  
}


bpwmissing <- function(t.est,present.data,n.com){
  for (h in 1:ncol(t.est)){
    if (length(which(t.est[,h]==0))!=0){
      t.est[,h][-which(t.est[,h]==0)] = present.data[t.est[,h]]
    } 
    else {t.est[,h] = present.data[t.est[,h]]}
  }
  diag(t.est) = n.com
  invisible(t.est)
}


