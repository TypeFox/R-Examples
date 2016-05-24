plot.dke.fun <-
function(x,
		
			main=NULL,
			sub = NULL, 
			xlab=NULL, 
			ylab=NULL,
			type="l",
                     	las=1,
			lwd=1,
			col="blue",
			lty=1,...)
                    {
    class(x) <- "dke.fun"
	
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)){
	                           ylab <- "density function"
	                }
    if(is.null(main)){ 
	     main <- "Kernel density estimate"
	                            
	                }
    if(is.null(sub)){
	     if(x$kernel=="GA") kerf="gamma"
	     	else if(x$kernel=="BE") 		kerf= "extended beta"
	        else if(x$kernel=="LN") 		kerf= "lognormal"
	        else if(x$kernel=="RIG") 		kerf= "reciprocal inverse Gaussian"

	         sub <- paste("Kernel = ",kerf,";", " h_n = ", formatC(x$h),";"," C_n = ",formatC(x$C_n))
	                }  

  
 	hist(x$data,xlab=xlab,ylab=ylab,sub=sub,probability=TRUE,main=main,border ="gray",... )
	lines(x$eval.points,x$est.fn,lwd=2,col=col,ylim=c(min(x$est.fn,na.rm = TRUE),max(x$est.fn,na.rm = TRUE)),...)

	
    invisible(NULL)
}
