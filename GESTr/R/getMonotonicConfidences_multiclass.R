#
# function for ensuring *sensible* monotonicity of multi-class models
#
###########################################################################


getMonotonicConfidences_multiclass <- function(ClassScoreList,exprs,verbose=FALSE){

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

	# correct low values that are too high

	newLoScores <- ClassScoreList[[1]]
	for(j in hiClass){
		if(sum(ClassScoreList[[j]]>0.1)>1){
			mixprops <- pnorm((exprs-exprs[which(ClassScoreList[[hiClass]]==max(ClassScoreList[[hiClass]]))][1])/sd(exprs[ClassScoreList[[j]]>0.1]),lower.tail=TRUE)
		}
		else{
			mixprops <- pnorm((exprs-exprs[which(ClassScoreList[[hiClass]]==max(ClassScoreList[[hiClass]]))][1])/sd(exprs),lower.tail=TRUE)
		}
		newLoScores <- newLoScores*(1-mixprops)
	}

	# correct high values that are too high
	newHiScores <- ClassScoreList[[hiClass]]
	for(j in 1){
                if(sum(ClassScoreList[[j]]>0.1)>1){
                        mixprops <- pnorm((exprs-exprs[which(ClassScoreList[[1]]==max(ClassScoreList[[1]]))][1])/sd(exprs[ClassScoreList[[j]]>0.1]),lower.tail=FALSE)
                }
                else{
                        mixprops <- pnorm((exprs-exprs[which(ClassScoreList[[1]]==max(ClassScoreList[[1]]))][1])/sd(exprs),lower.tail=FALSE)
                }
                newHiScores <- newHiScores*(1-mixprops)
        }

	# correct low values that are too low

	if(length(lo_indices)>0){
	if(verbose) cat("correcting 'lo' confidences \n")
		max_lo <- which(ClassScoreList[[1]]==max(ClassScoreList[[1]]))
		if(length(max_lo)>1){
			max_lo <- max_lo[1]
		}
		if(sum(ClassScoreList[[1]]>0.1)>1){
			mixprops <- pnorm((exprs-exprs[max_lo])/sd(exprs[ClassScoreList[[1]]>0.1]),lower.tail=FALSE)
		}
		else{
			mixprops <- pnorm((exprs-exprs[max_lo])/sd(exprs),lower.tail=FALSE)
		}
		for(j in 2:(hiClass)){
			newLoScores <- newLoScores + (mixprops*ClassScoreList[[j]])			
		}
		
	}
	newLoScores[newLoScores < 0] <- 0
	newLoScores[newLoScores > 1] <- 1
	
	# correct high values
	if(length(hi_indices)>0){
	if(verbose) cat("correcting 'hi' confidences \n")
		max_hi <- which(ClassScoreList[[hiClass]]==max(ClassScoreList[[hiClass]]))
		if(length(max_hi)>1){
			max_hi <- max_hi[length(max_hi)]
		}
		if(sum(ClassScoreList[[hiClass]]>0.1)>1){
			mixprops <- pnorm((exprs-exprs[max_hi])/sd(exprs[ClassScoreList[[hiClass]]>0.1]),lower.tail=TRUE)
		}
		else{
			mixprops <- pnorm((exprs-exprs[max_hi])/sd(exprs),lower.tail=TRUE)
		}
		for(j in 1:(hiClass-1)){
			newHiScores <- newHiScores + (mixprops*ClassScoreList[[j]])
		}
	}
	newHiScores[newHiScores < 0] <- 0
        newHiScores[newHiScores > 1] <- 1

	list(confidences_lo=newLoScores,confidences_hi=newHiScores)

}
