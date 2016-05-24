filter.absent <-
function(x,pct){
	present <- TRUE
	col <- length(x)/2
	pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
	if(pct.absent > pct){present <- FALSE}
	present
	}

