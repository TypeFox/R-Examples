#
# getMonotonicConfidences2: simple version of ensuring monotonic confidences
# (may be slower than other version)
#
#############################################################################

getMonotonicConfidences <- function(ClassScoreList){

        hiClass <- length(ClassScoreList)

        # find inversions for low-class
        lo_indices <- which((c(ClassScoreList[[1]],0)-c(1,ClassScoreList[[1]]))>0)
        if(ClassScoreList[[1]][1] < ClassScoreList[[1]][2]){
                lo_indices <- c(1,lo_indices)
        }
        if(ClassScoreList[[1]][length(ClassScoreList[[1]])] > ClassScoreList[[1]][length(ClassScoreList[[1]])-1]){
                lo_indices <- c(lo_indices,length(ClassScoreList[[1]]))
        }

        # find inversions for high-class
        hi_indices <- which((c(ClassScoreList[[hiClass]],1)-c(0,ClassScoreList[[hiClass]]))<0)
        if(ClassScoreList[[hiClass]][1] > ClassScoreList[[hiClass]][2]){
                hi_indices <- c(1,hi_indices)
        }
        if(ClassScoreList[[hiClass]][length(ClassScoreList[[hiClass]])] < ClassScoreList[[hiClass]][length(ClassScoreList[[hiClass]])-1]){
                hi_indices <- c(hi_indices,length(ClassScoreList[[hiClass]]))
        }

	# if inversions in low-class, ensure monotonicity point-wise

	newLoScores <- ClassScoreList[[1]]
	if(length(lo_indices)>0){
		pivot <- min(which(ClassScoreList[[1]]==min(ClassScoreList[[1]])))
		if(sum(lo_indices>pivot) > 0){
			# fix points to right of min value
			for(i in (pivot+1):length(ClassScoreList[[1]])){
				if(ClassScoreList[[1]][i] > ClassScoreList[[1]][i-1]){
					# try adding all descending scores from other classes at this point
					newLoScores[i] <- newLoScores[i-1]
					for(j in 2:hiClass){
						if(ClassScoreList[[j]][i] < ClassScoreList[[j]][i-1]){
							newLoScores[i] <- newLoScores[i] + (ClassScoreList[[j]][i] - ClassScoreList[[j]][i-1])
						}
					}
				}
				else{
					newLoScores[i] <- newLoScores[i-1] + (ClassScoreList[[1]][i]-ClassScoreList[[1]][i-1])
				}
				if(newLoScores[i] < 0){ newLoScores[i] <- 0 }
			}
		}
		if(sum(lo_indices<pivot) > 0){
			# fix points to left of min value
			for(i in (pivot-1):1){
				if(ClassScoreList[[1]][i] < ClassScoreList[[1]][i+1]){
					newLoScores[i] <- newLoScores[i+1]
					for(j in 2:hiClass){
						if(ClassScoreList[[j]][i] > ClassScoreList[[j]][i+1]){
							newLoScores[i] <- newLoScores[i] + (ClassScoreList[[j]][i] - ClassScoreList[[j]][i+1])
						}
					}
				}
				else{
					newLoScores[i] <- newLoScores[i+1] + (ClassScoreList[[1]][i] - ClassScoreList[[1]][i+1])
				}
				if(newLoScores[i] > 1){ newLoScores[i] <- 1 }
			}
		}
	}

	# if inversions in hi-class, ensure monotonicity point-wise

	newHiScores <- ClassScoreList[[hiClass]]
	if(length(hi_indices)>0){
		pivot <- max(which(ClassScoreList[[hiClass]]==max(ClassScoreList[[hiClass]])))
		if(sum(hi_indices>pivot) > 0){
			# fix points to right of max value
			for(i in (pivot+1):length(ClassScoreList[[hiClass]])){
				if(ClassScoreList[[hiClass]][i] < ClassScoreList[[hiClass]][i-1]){
					newHiScores[i] <- newHiScores[i-1]
					for(j in 1:(hiClass-1)){
						if(ClassScoreList[[j]][i] > ClassScoreList[[j]][i-1]){
							newHiScores[i] <- newHiScores[i] + (ClassScoreList[[j]][i] - ClassScoreList[[j]][i-1])
						}
					}
				}
				else{
					newHiScores[i] <- newHiScores[i-1] + (ClassScoreList[[hiClass]][i] - ClassScoreList[[hiClass]][i-1])
				}
				if(newHiScores[i] > 1){ newHiScores[i] <- 1 }
			}
		}
		if(sum(hi_indices<pivot) > 0){
			# fix points to left of min value
			for(i in (pivot-1):1){
				if(ClassScoreList[[hiClass]][i] > ClassScoreList[[hiClass]][i+1]){
					newHiScores[i] <- newHiScores[i+1]
					for(j in 1:(hiClass-1)){
						if(ClassScoreList[[j]][i] < ClassScoreList[[j]][i+1]){
							newHiScores[i] <- newHiScores[i] + (ClassScoreList[[j]][i] - ClassScoreList[[j]][i+1])
						}
					}
				}
				else{
					newHiScores[i] <- newHiScores[i+1] + (ClassScoreList[[hiClass]][i] - ClassScoreList[[hiClass]][i+1])
				}
				if(newHiScores[i] < 0){ newHiScores[i] <- 0 }
			}
		}
	}

	list(confidences_lo=newLoScores,confidences_hi=newHiScores)

}

