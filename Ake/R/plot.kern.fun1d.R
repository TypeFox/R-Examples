plot.kern.fun1d <-
function(x,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      pch=18,las=1,lwd=1,...)
                  {
    class(x) <- "kern.fun"
    kernel <- x$ker
    if(is.null(xlab)) xlab <- "y"
    if(is.null(ylab)) ylab <- "Prob(y)" 
    if(is.null(main)){ 
    	    
    	                 main <- "Kernel function"
    	                }               
    
    plot.default(x$t,x$kx,pch=pch,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)

}
