Vplot <- function(object, units=NULL, ...)  {
    nunits <- nrow(object$V)
    if(is.null(units)) {
        ### plot 5 Valpha curves
        units <- sample(c(1:nunits),size=5)
    }
    nplotunits <- length(units)
    
    par(mar=c(5.1,4.6,1.6,2.1))
    plot(object$aux$alpha.grid,object$aux$Vmarginals,type="n",xlim=c(1/nunits,1 - 1/nunits),ylim=c(0,1),
         xlab = expression(alpha),ylab=expression(V[alpha]),cex.lab=1.3)
    lines(object$aux$alpha.grid,object$aux$Vmarginals,lwd=2,lty=2)
    for(k in 1:nplotunits) {
        lines(object$aux$alpha.grid,object$aux$V[units[k],], lwd=2)
    }
    legend("topleft", legend=c(expression(lambda(alpha)),
           expression(V[alpha](D[i])) ),lwd=2, lty=c(1,2), bty="n",cex=1.4)
} 
