get.cycle.durations <-
function(dataset){
	ncycles <- dim(dataset)[2]
	durations <- NULL
	for(i in 1:ncycles){
		set1 <- rm.end.NAs(dataset[,i])
		durations[i] <- length(set1)
		}
	return(durations)
	}

