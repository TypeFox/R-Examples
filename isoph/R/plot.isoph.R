plot.isoph=function(x, which=1, lty=1, lcol=1, lwd=1, pch=19, pcol=1, pcex=0.7, main=NULL, ylab=NULL, xlab=NULL, lglab=NULL, lgloc=NULL, lgcex=0.9, ylim=NULL, xlim=NULL, ...){
  y.obs=x$psi
  z.obs=x$z
  z.range=x$z.range
  hr.obs=exp(y.obs)
  
  n=length(y.obs)
  y=c(y.obs[1],y.obs,y.obs[n])
  hr=c(hr.obs[1],hr.obs,hr.obs[n])
  z=c(z.range[1],z.obs,z.range[2])
  
  if(is.null(main)) main=paste("Isotonic proportional hazards model\n(monotone ",x$shape," covariate effect)",sep="")
  if(is.null(xlab)) xlab=x$formula[[3]]
  if(is.null(xlim)) xlim=z.range
  
  if(is.null(lglab)) lglab="Potential jump point"
  if(is.null(lgloc)){
    if(x$shape=='increasing'){       lgloc='topleft'
    }else if(x$shape=='decreasing'){ lgloc='topright' 
    }
  }  
  
  if(which==1){
    if(is.null(ylab)) ylab=expression(hat(psi))
    if(is.null(ylim)) ylim=range(y)
    
    plot(y~z,type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=lcol, lty=lty, lwd=lwd)
    if( is.character(pch) || (is.numeric(pch) & is.finite(pch)) )
      points(y.obs~z.obs, pch=pch, col=pcol, cex=pcex)
  }else if(which==2){  
    if(is.null(ylab)) ylab=expression(exp(hat(psi)))
    if(is.null(ylim)) ylim=range(exp(y))    
    
    plot(hr~z,type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=lcol, lty=lty, lwd=lwd)
    if( is.character(pch) || (is.numeric(pch) & is.finite(pch)) )
      points(hr.obs~z.obs, pch=pch, col=pcol, cex=pcex)
  }
  
  #legends
  if(!is.na(lglab))
    legend(x=lgloc, legend=lglab, pch=pch, col=pcol, bty='n', cex=lgcex)  
}