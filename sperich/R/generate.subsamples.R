generate.subsamples <-
function(number.of.occurrences, fold, loocv.limit){

	if (number.of.occurrences < loocv.limit){
		#loocv
		subsamples <- matrix(0,number.of.occurrences, number.of.occurrences-1)
		for (i in 1:number.of.occurrences){
			if ((i > 1)&&(i < number.of.occurrences)){
				subsamples[i,] <- c(1:(i-1),(i+1):number.of.occurrences)
			}
			if (i==1){
				subsamples[i,] <- (i+1):number.of.occurrences
			}
			if (i == number.of.occurrences){
				subsamples[i,] <- 1:(number.of.occurrences-1)
			}
		}
	} else {
		#x-fold cross validation
		members <- 1:number.of.occurrences
		number.of.members <- ceiling(number.of.occurrences/fold)
		groups <- matrix(0, fold, number.of.members)

		for (i in 1:fold){
			for (j in 1:number.of.members){
				if (length(members)==0){
					break
				}
				chosen <- round(runif(1,0,length(members)))
				if (chosen == 0){
					chosen <- length(members)
				}
				groups[i,j] <- members[chosen]
				members <- members[which(members != members[chosen])]
			}
			groups[i,] <- sort(groups[i,])
		}

		if (length(which(groups[fold,]!=0))==0){
			groups <- groups[1:(fold-1),]
		}
		subsamples <- matrix(0, dim(groups)[1], (dim(groups)[1]-1)*number.of.members)

		for (i in 1:dim(groups)[1]){
			subsample <- c()
			for (j in 1:dim(groups)[1]){
				if (j != i){
					subsample <- c(subsample, groups[j,])
				}
			}
			subsamples[i,] <- sort(subsample)
		}
	}

	return(subsamples)
}