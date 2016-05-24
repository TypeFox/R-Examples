mway.fact <- function(pnom,plab,perc.v,levels){

 npar <- length(pnom)
 pv <- matrix(nrow=levels,ncol=npar); pf <- list(); pfr <- list()
 n <- which(seq(3,17,2)==levels)
 for(i in 1:npar){ 
  v <- perc.v*pnom[i]/100
  pv[,i] <- seq(pnom[i]-n*v, pnom[i]+n*v,v)
  pf[[i]] <- factor(pv[,i])
 } 
 pval <- matrix(nrow=levels^npar,ncol=npar)
 for(i in 1:npar){ 
  pfr[[i]] <- gl(levels,levels^(i-1),levels^npar,labels=pv[,i])
  pval[,i] <- pv[c(pfr[[i]]),i]
 }
 return(list(plab=plab, pval=pval, pnom=pnom, fact=T, pv=pv, pf=pf, pfr=pfr))
}

