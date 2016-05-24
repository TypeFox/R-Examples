plot.sparsenet=function(x, xvar=c("rsq","lambda","norm"),which.gamma=NULL,label=FALSE,...){
  oldpar=par(mar=c(4,4,4,1))
  on.exit(par(oldpar))
  xvar=match.arg(xvar)
 switch(xvar,
        rsq={iname="Fraction Training Variance Explained";xvar="dev"},
        lambda={iname=expression(log (lambda))},
        norm={iname="L1 Norm"}
        )
  lamax=x$max.lambda
  coeflist=x$coefficients
  ngamma=length(coeflist)
  if(is.null(which.gamma))which.gamma=1:ngamma
  coeflistseq=seq(along=coeflist)
  which.gamma=coeflistseq[match(which.gamma,coeflistseq,0)]
  rsq=x$rsq
  ylims=range(sapply(coeflist[which.gamma],function(x)range(x$beta)))
  for(i in which.gamma){
      x=coeflist[[i]]
      beta=x$beta
      p=nrow(beta)
      beta=cbind2(rep(0,p),beta)
      plotCoef(beta,lambda=c(lamax,x$lambda),df=c(0,x$df),dev=c(0,rsq[i,]),label=label,xvar=xvar,xlab=iname,ylim=ylims,...)
      
      if(x$gamma> 1000)mlab="Lasso"
      else if (x$gamma<1.001)mlab="Subset"
      else mlab=paste("Gamma =",format(round(x$gamma,1)))
mtext(mlab,3,2)
  }
invisible()
}
