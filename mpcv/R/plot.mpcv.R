plot.mpcv <- function(x, ...)
{
  n.integr=100
  n.features <- ncol(x$x)
  ForIds <- c(1:n.features)
  ForIds <- ForIds[-x$indepvar]
  
  marginal.median <- apply(x$x, 2, median)
  
  points.integr <- seq(min(x$x[,x$indepvar]), max(x$x[,x$indepvar]), length=n.integr)
  models.integr.u <- matrix(ncol=n.features, nrow=length(points.integr))
  models.integr.l <- matrix(ncol=n.features, nrow=length(points.integr))
  
  a.u <- x$coef.u
  a.l <- x$coef.l
  
  for(i in ForIds)
  {
    models.integr.u[,i] <- a.u[i,1]-a.u[i,2]+a.u[i,3]*points.integr-a.u[i,4]*points.integr-a.u[i,5]*points.integr^2
    models.integr.l[,i] <- a.l[i,1]-a.l[i,2]+a.l[i,3]*points.integr-a.l[i,4]*points.integr+a.l[i,5]*points.integr^2
  }
  
  
  
  graph.max <- matrix(ncol=n.features)
  graph.min <- matrix(ncol=n.features)
  graph.max[x$indepvar] <- x$USL[x$indepvar]
  graph.min[x$indepvar] <- x$LSL[x$indepvar]
  
  par(mfrow=c(1,n.features-1), oma=c(0,0,2,0), xpd=NA)
  # values of models of each variable and a plot
  for(i in ForIds)
  {
    graph.max[i] <- max(models.integr.u[,i], x$USL[i])
    graph.min[i] <- min(models.integr.l[,i], x$LSL[i])
    par(mar=c(5.1, 5.1, 4.1, 2.1)) 
    plot(x$x[,x$indepvar],x$x[,i], type="p", pch=21, cex=1, col="black",  bg="darkgrey", xlim=c(.95*x$LSL[x$indepvar],1.05*x$USL[x$indepvar]), 
         ylim=c(.95*graph.min[i], 1.05*graph.max[i]),
         xlab=colnames(x$x)[x$indepvar], ylab=colnames(x$x)[i],
         cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1, mar=c(5,10,2,2))
    rect(x$LSL[x$indepvar], x$LSL[i], x$USL[x$indepvar], x$USL[i], col="darkgreen", density=0, lwd=2, lty = 2)
    lines(points.integr,models.integr.u[,i], col="red", lwd=2)
    lines(points.integr,models.integr.l[,i], col="red", lwd=2)
    lines(c(points.integr[1], points.integr[1]), c(models.integr.l[1,i], models.integr.u[1,i]), , col="red", lwd=2)
    lines(c(points.integr[n.integr], points.integr[n.integr]), c(models.integr.l[n.integr,i], models.integr.u[n.integr,i]), , col="red", lwd=2)
    points(x$Target[x$indepvar],x$Target[i], col="darkgreen",pch=10, cex=3, lwd=2)
    points(marginal.median[x$indepvar],marginal.median[i], col="red",pch="+", cex=2)
  }
  legend("topright", inset=-.1, ncol=3, legend=c("Specification limits", "Target", "Median"), bty="n",
         col=c("darkgreen", "darkgreen", "red"), lty = c(2,NA,NA), lwd=2, pch = c(NA, 10, 43), pt.cex=2, cex=1, merge=FALSE) 
  
}
