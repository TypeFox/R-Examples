defineAxis<-function(major.ticks, limits=NULL, minor.ticks=NULL, major.tick.labels=major.ticks, both.sides=FALSE){
	l <- list()
	
	l$limits <- NULL
	if (!is.null(limits)){
		l$limits <- limits
	}
	
	l$major.ticks <- major.ticks
	l$major.tick.labels <- major.tick.labels
	l$both.sides <- both.sides
	
	#remove any duplicates of primary ticks in seconday ticks
	if(!is.null(minor.ticks)){
		cull <- NULL
		for(i in 1:length(minor.ticks)){
			if(any(abs(minor.ticks[i]-major.ticks) < .Machine$double.eps ^ 0.5))
				cull<-cbind(i,cull)			
		}
		if(!is.null(cull))l$minor.ticks <- minor.ticks[-cull]
	}
	else
		l$minor.ticks <- NULL

    return(l)

}
