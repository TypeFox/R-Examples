get.breaks.vector <-
function(CycleBreaksMatrix){
	breaks.vector <- NULL
	for(i in 1:dim(CycleBreaksMatrix)[2]){
		set <- na.omit(CycleBreaksMatrix[,i])
		lset <- length(set)-1
		breaks.vector <- append(breaks.vector, set[1:lset])
		}	
	return(breaks.vector)
	}

