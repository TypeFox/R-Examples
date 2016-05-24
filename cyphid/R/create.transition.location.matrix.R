create.transition.location.matrix <-
function(transitionVector, cycleCounts){
	MaxTime <- max(cycleCounts)
	end <- 0
	transMat <- NULL
	for(i in 1:length(cycleCounts)){
		start <- end + 1
		end <- start + cycleCounts[i] - 1
		set1 <- transitionVector[start:end]
		transMat <- add.col(transMat, set1, MaxTime)
		}
	return(transMat)
	}

