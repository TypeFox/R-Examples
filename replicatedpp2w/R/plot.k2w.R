plot.k2w<- function (x, trat=NULL, ...,  lty = NULL, col = NULL, 
    lwd = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
     legend = TRUE, legendpos = "topleft",  fun="L", main=NULL) 
{
        if(is.null(main)){
	   main <- if (is.language(substitute(x))) 
             short.deparse(substitute(x))
           else ""
        }
        if (is.null(lwd)) lwd <- c(1,1)
	if(is.null(xlab)) xlab <- "r"


       if(is.null(trat)){
          if(is.null(x$tratB)){
	        if (is.null(col)) col <- unique(unclass(x$tratA))+1
		if(is.null(lty)) lty <- unique(unclass(x$tratA))
		leveltrat <- levels(x$tratA)
		titulo <- x$nameA
		
		
		if(fun=="L") {
		if (is.null (ylab)) ylab <- "L(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, sqrt(x$KrepA/pi)-x$r, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab,xlim=xlim, ylim=ylim,  main=main, ... )
		   lines(x$r, sqrt(x$K0i/pi)-x$r, lwd=lwd[2])
	       }
	       if(fun!="L") {
		  if (is.null (ylab)) ylab <- "K(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, x$KrepA, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, x$K0i, lwd=lwd[2])
	       }
	       
	       
	 }
	   
	 if(!is.null(x$tratB)){
		 if (is.null(col)) col <- unique(unclass(x$tratAB))+1
		if(is.null(lty)) lty <- unique(unclass(x$tratAB))
		leveltrat <- levels(x$tratAB)
		titulo <- "Interaction"
		
		if(fun=="L") {
		if (is.null (ylab)) ylab <- "L(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, sqrt(x$KrepAB/pi)-x$r, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, sqrt(x$K0ij/pi)-x$r, lwd=lwd[2])
	       }
	       if(fun!="L") {
		  if (is.null (ylab)) ylab <- "K(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, x$KrepAB, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, x$K0ij, lwd=lwd[2])
	       }
	 }
      }
      if(!is.null(trat)){ 
        if(trat =="tratA"){
	        if (is.null(col)) col <- unique(unclass(x$tratA))+1
		if(is.null(lty)) lty <- unique(unclass(x$tratA))
		leveltrat <- levels(x$tratA)
		titulo <- x$nameA
		
		if(fun=="L") {
		if (is.null (ylab)) ylab <- "L(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, sqrt(x$KrepA/pi)-x$r, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, sqrt(x$K0i/pi)-x$r, lwd=lwd[2])
	       }
	       if(fun!="L") {
		  if (is.null (ylab)) ylab <- "K(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, x$KrepA, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, x$K0i, lwd=lwd[2])
	       }
	       
	       
	  }
      
      
        if(trat =="tratB"){
	        if(is.null(x$tratB)) stop (paste("only tratA , i.e.,", x$nameA,"employed in this analysis","\n"))
	        if (is.null(col)) col <- unique(unclass(x$tratB))+1
		if(is.null(lty)) lty <- unique(unclass(x$tratB))
		leveltrat <- levels(x$tratB)
		titulo <- x$nameB
		
		if(fun=="L") {
		if (is.null (ylab)) ylab <- "L(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, sqrt(x$KrepB/pi)-x$r, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, sqrt(x$K0j/pi)-x$r, lwd=lwd[2])
	       }
	       if(fun!="L") {
		  if (is.null (ylab)) ylab <- "K(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, x$KrepB, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, x$K0j, lwd=lwd[2])
	       }
	       
	       
	 }
      
        if(trat=="tratAB"){
	        if(is.null(x$tratAB)) stop (paste("only tratA , i.e.,", x$nameA,"employed in this analysis","\n"))
		if (is.null(col)) col <- unique(unclass(x$tratAB))+1
		if(is.null(lty)) lty <- unique(unclass(x$tratAB))
		leveltrat <- levels(x$tratAB)
		titulo <- "interaction"
		
		if(fun=="L") {
		if (is.null (ylab)) ylab <- "L(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, sqrt(x$KrepAB/pi)-x$r, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, sqrt(x$K0ij/pi)-x$r, lwd=lwd[2])
	       }
	       if(fun!="L") {
		  if (is.null (ylab)) ylab <- "K(r)"
		   #matplot(x$r, sqrt(x$dataKijk/pi)-x$r, type="l", col=unclass(x$tratA)+1, lty=unclass(x$tratA))
		   matplot(x$r, x$KrepAB, type="l", col=col, lty=lty, lwd=lwd[1], xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, ... )
		   lines(x$r, x$K0ij, lwd=lwd[2])
	       }
	 }
      }
      
        if(legend ){
	   legend(legendpos, legend=c(leveltrat, "Global"), 
	                col=c(col,1), lty=c(lty,1), title = titulo,
			lwd=c(rep(lwd[1], length(leveltrat)), lwd[2]))
       }

}
      