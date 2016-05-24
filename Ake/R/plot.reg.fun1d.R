plot.reg.fun1d <-
function(f,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      pch=18,las=1,lwd=1,...)
                  {
    class(f) <- "reg.fun"
    kernel <- f$ker
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- "y" 
    if(is.null(main)){ 
    	    
    	                 main <- "kernel regression"
    	                }               
    
    plot.default(f$data,f$y,pch=pch,col="grey",las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)


}
