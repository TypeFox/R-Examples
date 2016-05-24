plot.kpmfe.fun1d <-
function(x,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="h",las=1,lwd=1,...)
                  {
    class(x) <- "kpmfe.fun"
    kernel <- x$ker
    if(is.null(xlab)) xlab <- "y"
    if(is.null(ylab)) ylab <- "Prob(y)" 
    if(is.null(main)){ 
    	    
    	                 main <- "Discrete kernel estimaton"
    	                }               
    
    plot.default(x$eval.points,x$est.fn, type="h",las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)

}
