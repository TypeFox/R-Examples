get.cycle.counts <-
function(closeMatrix){
	cycle.counts <- NULL
	for(i in 1:dim(closeMatrix)[2]){
		ncycles <- length(na.omit(closeMatrix[,i]))
		cycle.counts <- append(cycle.counts, ncycles)
		}
	return(cycle.counts)
	}

