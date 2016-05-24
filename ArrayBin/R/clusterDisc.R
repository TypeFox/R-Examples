clusterDisc <- function(x,use.gap){

	outvals <- rep(NA,length(x))
	notNAs <- which(!is.na(x))
	if(length(notNAs)>1){
		x <- x[notNAs]
		clusters <- list()
		gap.stats <- list()
		if(use.gap){
			for (k in 1:2){
				clusters[[k]] <- mskmeans(x,k)
				gap.stats[[k]] <- gap(data=as.matrix(x),class=clusters[[k]]$cluster,B=100,cluster.func=mskmeans)
			}
			# use gap statistic to choose k=1 or k=2
			k <- ifelse(gap.stats[[1]][1]>=(gap.stats[[2]][1]-gap.stats[[2]][2]),1,2)
		}
		if(!use.gap){
			k <- 2
			clusters[[k]] <- mskmeans(x,k)
		}
		outvals[notNAs] <- clusters[[k]]$cluster-1
	}
	outvals

}
