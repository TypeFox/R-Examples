unisegmixt <- function(Y,CGHo,Kmax,phi){
  P            = CGHo["nblevels"]
  n.com        = length(Y)
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  out          = ClassiSeg(x,phi[1:P],Kmax)
	
  loglik = sapply(1:Kmax,FUN=function(k){
    th        = out$t.est[k,1:k]
    rupt      = matrix(ncol=2,c(c(1,th[1:k-1]+1),th))    
    return(lvinc(Y,phi,rupt,P))
  })
  
  t.est     = bpwmissing(out$t.est,present.data,n.com)
  th        = t.est[Kmax,1:Kmax]
  rupt      = matrix(ncol=2,c(c(1,th[1:Kmax-1]+1),th))
  mu        = data.frame(begin = rupt[,1],
    end   = rupt[,2],
    mean  = rep(NA,Kmax)
    )  
	
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


lvinc  <- function (Y,phi,rupt,P){
  x           =  Y
  logdensity  = t(apply(rupt, 1,FUN = function(y){
    xk  = x[y[1]:y[2]]
    xk  = xk[!is.na(xk)]
    invisible(logdens(xk,P, phi))
  }))
  K       = nrow(logdensity)
  P       = ncol(logdensity)
  tau     = sapply(1:P,FUN = function(p){logdensity[,p]+log(phi[p+2*P])})
  tau     = matrix(tau,ncol=P)
  tau_max = apply(tau,1,max)
  tau     = exp(tau-matrix(rep(tau_max,P),ncol=P))
  lvinc   = sum(log( apply(tau,1,sum)) + tau_max)
  return(lvinc)
}
