getPeaks <-
function(x,y, sep=8,thr=0.01){
	x <- c(0,0,0,0,0,x,0,0,0,0,0)
	y <- c(0,0,0,0,0,y,0,0,0,0,0)
	peak.num <- 0
	num <- -1
	peak.max <- c()
	peak.vmax <- c()
	x.length <- length(x)
	while(1){
		x.step <- as.integer(x.length/sep)
		num <- peak.num
		peak.values <- c()
		peak.xpos <- c()
		i <- 1
		sep <- sep*2
		if(sep >= length(x))break
		while(i <= x.length-x.step-1){
			loc.max <- max(y[i:(i+x.step)])
			if(abs(loc.max)<=thr){
				i <- i+x.step			
				next}
			loc.x <- x[i:(i+x.step)]
			max.x <- loc.x[y[i:(i+x.step)]==loc.max]
			if((loc.max != y[i])&(loc.max != y[i+x.step])){
				peak.values <- cbind(peak.values,loc.max)
				peak.xpos <- cbind(peak.xpos,max.x)
			}
			if((loc.max == y[i])&(i>=5)){
				if(loc.max==max(y[(i-min(x.step/2,4)):(i+min(x.step/2,4))])){
					peak.values <- cbind(peak.values,loc.max)
					peak.xpos <- cbind(peak.xpos,max.x)
				}
			}
			if((loc.max == y[i+x.step])&(i<=(length(x)-5))){
				if(loc.max==max(y[(i+x.step-min(x.step/2,4)):(i+x.step+min(x.step/2,4))])){
					peak.values <- cbind(peak.values,loc.max)
					peak.xpos <- cbind(peak.xpos,max.x)
				}
			}
			i <- i+x.step+1
		}
		peak.xpos <- as.character(peak.xpos)
		peak.values <- as.character(peak.values)
		if(peak.num <= length(unique(peak.xpos))){
			peak.num <- length(unique(peak.xpos))
			peak.max <- as.numeric(peak.xpos)
			peak.vmax <- as.numeric(peak.values)
		}
	}	
	return(cbind(peak.max,peak.vmax))
}
