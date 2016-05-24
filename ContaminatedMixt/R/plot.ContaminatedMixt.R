plot.ContaminatedMixt <- function(x, criterion="BIC",contours=FALSE, xmarg=1, ymarg=2, res=200, levels=seq(.0001,1,by=0.01), ...){


  criterion <- match.arg(criterion,.ICnames())
  bivres <- getBestModel(x,criterion=criterion)$models[[1]]
  
  
  plot(bivres$X[,c(xmarg,ymarg)], col="white", ...)
  if(contours){
    lims <- par()$usr
    xseq <- seq(lims[1], lims[2], length.out=res)
    yseq <- seq(lims[3], lims[4], length.out=res)
    resgood <- resbad <- rescont <- array(0,c(res,res,bivres$G))
    val     <- array(0,c(res,res))
    for(g in 1:bivres$G){
      resgood[,,g] <- outer(xseq,yseq,.bivnorm,bivres$mu[c(xmarg,ymarg),g],bivres$Sigma[c(xmarg,ymarg),c(xmarg,ymarg),g])
      resbad[,,g]  <- outer(xseq,yseq,.bivnorm,bivres$mu[c(xmarg,ymarg),g],bivres$eta[g]*bivres$Sigma[c(xmarg,ymarg),c(xmarg,ymarg),g])
      rescont[,,g] <- bivres$alpha[g]*resgood[,,g]+(1-bivres$alpha[g])*resbad[,,g]
      val          <- val + bivres$prior[g]*rescont[,,g]
    }
    contour(x=xseq, y=yseq, z=val, add=TRUE, levels=levels, col=rgb(0.5,0.5,0.5,alpha=0.7))
  }
  text(bivres$X[,c(xmarg,ymarg)], labels=bivres$detection$innergroup, col=bivres$group, ...)
  
}
