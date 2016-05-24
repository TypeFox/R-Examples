kts.pair <- function (dat, grp, k, display = TRUE, length = 40, med = FALSE){ 

	labels <- as.character(unique(grp)[order(unique(grp))])

	if(!is.null(k)){
		if(k>length){stop("The value of length has to be at least as big as k")}
		if(k%%2==0){stop("k needs to be an odd number")}
	}


	if(med){
		median_0 <- c(apply(dat[,which(grp==1)], 1, function(x) median(x,na.rm=TRUE)))
		median_1 <- c(apply(dat[,which(grp==0)], 1, function(x) median(x,na.rm=TRUE)))
		mean_med <- (median_0+median_1)/2
		dat <- dat-mean_med
	}

	replacena <- floor(min(dat[is.na(dat)==FALSE]))-100000
		
	## Call the C code to compute the k best pairs.
	out <- .Call("cktspair", as.double(rank_na(dat,replacena)), as.double(grp),as.integer(k), as.double(replacena), as.integer(length), PACKAGE="ktspair")

	index <- matrix(ncol=2, nrow=k)
	rankscore <- c()
	ktspscore <- c()
		
	for (i in 1:k){
		if((out[[1]][3 + (i-1) * 4]!=0) && (out[[1]][i * 4]!=0)){
			index[i,1] <- out[[1]][3 + (i-1) * 4]
			index[i,2] <- out[[1]][i * 4]
			rankscore[i] <- out[[1]][2 + (i-1) * 4]
			ktspscore[i] <- out[[1]][1 + (i-1) * 4]
		}
	}
	
	if(is.na(index[k,1])==TRUE){## Check if it was possible to compute k pairs, if not update the value of k.
		k2 <- min(which(is.na(index[,1])))-1
		if(display){
			cat("The k-TSP cannot be computed for k = ", k, "\n")
		}
		if(k2%%2==0){k<-k2-1}
		if(k2%%2==1){k <- k2}		
		if(display){
			cat("The value of k has been reduced to k = ", k, "\n")
		}
	}
	index <- matrix(c(index[1:k,1],index[1:k,2]),ncol=2, nrow=k) 

	ktspdat <- dat[c(index[,1],index[,2]),]
	accuracy_k <- NULL
	accuracy <- NULL
	sensitivity <- NULL
	specificity <- NULL

	if(med) med <- mean_med[c(index)]

	## Create a k-TSP object.
	ktsp <- list(index = index, ktspscore = ktspscore[1:k], grp = grp, ktspdat = ktspdat, k = k, labels =labels, rankscore = rankscore[1:k], accuracy =accuracy, accuracy_k = accuracy_k, sensitivity = sensitivity, specificity = specificity, med = med)
	class(ktsp) <- "ktsp"
	return(ktsp)
}

