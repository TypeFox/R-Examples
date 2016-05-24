plot.hcvd.fun1d <-
function(f,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                  {
    class(f) <- "CV.fun"
    kernel <- f$ker
    if(is.null(xlab)) xlab <- "h"
    if(is.null(ylab)) ylab <- "CV(h)" 
    if(is.null(main)){ 
    	    
    	                 main <- "Cross-validation"
    	                }               
    
    plot.default(f$seq_h,f$CV, type="l",las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)

}
