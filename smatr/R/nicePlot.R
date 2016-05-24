nicePlot<-function(xaxis,
                   yaxis,
                   log='', 
                   ann=par("ann"), 
                   xlab=NULL, 
                   ylab=NULL, 
                   tck=par("tck"),
                   frame.plot = TRUE, 
                   ... ){
  
	makeAxis<-function(axis, number, print.labels=TRUE, tck=0.25, ...){
		if(print.labels)
			axis(number, at=axis$major.ticks, labels=axis$major.tick.labels, tck=tck, ...)
		else
			axis(number, at=axis$major.ticks, labels=FALSE, tck=tck,...)
		if(length(axis$minor.ticks)>0)
    		axis(number, at=axis$minor.ticks, labels=FALSE,tck=tck*0.5, ...)
	}

	
	#check sufficient arguments present
	if(is.null(yaxis$limits) | is.null(yaxis$major.ticks)){
		cat("Data missing for Yaxis, check limits and major ticks\n\n");
	} else if(is.null(xaxis$limits) | is.null(xaxis$major.ticks)){
		cat("Data missing for Xaxis, check limits and major ticks\n\n");
	} else{
		
		#make empty plot
		plot(1, 1, type="n", log=log, axes=FALSE,ann=FALSE, ylim=yaxis$limits, 
			xlim=xaxis$limits, xaxs="i", yaxs="i")
	
		#plot axes
		makeAxis(yaxis, 2, 1, tck=tck, ...)
		makeAxis(xaxis, 1, 1, tck=tck, ...)
		
		#add axes on opposite side if appropriate
		if(yaxis$both.sides==TRUE)
			makeAxis(yaxis, 4, 0, tck=tck,...)	
		if(xaxis$both.sides==TRUE)	
			makeAxis(xaxis, 3, 0, tck=tck,...)
		#add box if appropriate
		if(frame.plot) 
		    box(...)
	
		#plot axis titles if appropriate
		if(ann && !is.null(xlab))
			mtext(xlab, side = 1, line = par("mgp")[1], outer = FALSE, at= NA)
		if(ann && !is.null(ylab))
			mtext(ylab, side = 2, line =  par("mgp")[1], outer = FALSE, at= NA)
	}
}
