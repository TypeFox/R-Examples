`plot.angledist` <-
function(x, ylim=NULL, add=FALSE, linecol="blue",xlab=expression(Leaf~angle~~(""^"o")),ylab="Density",main=NA,...){

	object <- x
	a <- object$fitdata
	if(is.na(main))main <- object$distribution
	
	if(!all(is.na(a))){
	
		if(max(a) < pi/2 + 0.01){
			warning("Converting from radians to degrees\n")
			a <- a * 180/pi
		}
	
		if(!add){hist(a, breaks=seq(0,90, by=5),
			freq=FALSE, 
			main=main,
			xlab=xlab, 
			ylab=ylab, ylim=ylim ,...)
		}
	
		curve(ftheta(x, distribution=object$distribution, distpars=object$distpars, degrees=TRUE), 
			add=TRUE, col=linecol)
			
		box()
	
	} else {   # No fit data, distribution was therefore defined with 'angledist' function.
	
		curve(ftheta(x, distribution=object$distribution, distpars=object$distpars, degrees=TRUE), 
			main=main, 
			xlab=xlab, 
			ylab=ylab, ylim=ylim,add=add,
			from=0, to=90,
			col=linecol)
		box()
	}

}

