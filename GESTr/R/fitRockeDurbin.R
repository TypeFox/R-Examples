#
# function to return list of parameters for Rocke-Durbin error model
# *** remember to use correct (log,exp) transform of raw values ***
# returns list of alpha, sd_epsilon and sd_eta
#
# ('theta' is similarity threshold for calling replicate groups)
##########################################################################

fitRockeDurbin <- function(x,theta){

	distMatrix <- as.matrix(dist(t(x)))
	cat("calculating Rocke-Durbin error model parameters \n")
	cat("finding replicate groups \n")
	replicateGroups <- findReplicateGroups(distMatrix,theta)
	cat(paste(length(replicateGroups),"groups in total \n"))
	cat("finding 'background' genes for each replicate group \n")

	findLowGenes <- function(samples,x){
	
		lowgenes <- sort(rowMeans(x[,samples]),decreasing=FALSE,index.return=TRUE)$ix[1:(dim(x)[1]/10)]
		newGenes <- lowgenes
		while(length(newGenes)>0 & length(lowgenes)>100){
			sample_mean <- median(x[lowgenes,samples])
			sample_sd <- mad(x[lowgenes,samples])/0.6745
			newLowGenes <- lowgenes[apply(x[lowgenes,samples],MARGIN=1,max)<(sample_mean+(2*sample_sd))]
			print(paste(length(newLowGenes),"still considered low"))
			newGenes <- setdiff(lowgenes,newLowGenes)
			lowgenes <- newLowGenes
		}
		cat(paste("converged onto set of",length(lowgenes),"genes \n"))
		lowgenes
	}

	lowgenes_list <- lapply(replicateGroups,findLowGenes,x)
	
	r <- length(replicateGroups)
	ljs <- unlist(lapply(replicateGroups,length))
	replicateGroups <- replicateGroups[ljs>1]
	ljs <- unlist(lapply(replicateGroups,length))
	r <- length(replicateGroups)
	l <- sum(ljs)

	cat("calculating parameters for low expression \n")

	alphas <- rep(NA,r)
	sd_epsilons <- rep(NA,r)
	for(i in 1:r){
		if(length(lowgenes_list[[i]])>1){
			relevantData <- 2^x[lowgenes_list[[i]],replicateGroups[[i]]]
			sds <- apply(relevantData,MARGIN=1,sd)
			means <- rowMeans(relevantData)
			alphas[i] <- mean(means)
			sd_epsilons[i] <- sqrt((1/length(sds))*sum(sds^2))
		}
	}
	
	# summarise alpha, sd_epsilon over all replicate groups
	# average weighted by number of samples in replicate group
	toRemove1 <- unique(c(which(is.na(alphas)),which(is.na(sd_epsilons))))
	if(length(toRemove1)>0){
		alpha <- (1/(l-r))*sum(alphas[-toRemove1]*(ljs[-toRemove1]-1))
	}
	if(length(toRemove1)==0){
		alpha <- (1/(l-r))*sum(alphas*(ljs-1))
	}
	cat(paste("estimate for background alpha=",alpha,"; log2scale alpha=",log(alpha,base=2)),"\n")
	if(length(toRemove1)>0){
		sd_epsilon <- sqrt((1/(l-r))*sum((sd_epsilons[-toRemove1]^2)*(ljs[-toRemove1]-1)))
	}
	if(length(toRemove1)==0){
		sd_epsilon <- sqrt((1/(l-r))*sum((sd_epsilons^2)*(ljs-1)))
	}

	findHighGenes <- function(samples,x,prop=0.1){
		order(rowMeans(x[,samples]),decreasing=TRUE)[1:round(nrow(x)*prop)]
	}

	cat("finding 'high-level' expressed genes \n")
	highgenes_list <- lapply(replicateGroups,findHighGenes,x=x)

	cat("calculating parameters for high expression \n")
	sd_etas <- rep(NA,r)
	for(i in 1:r){
		if(length(highgenes_list[[i]])>1){
		sds <- apply(log((2^x[highgenes_list[[i]],replicateGroups[[i]]])-alpha,base=2),MARGIN=1,sd)
		sds <- sds[!is.na(sds)]
		sd_etas[i] <- sqrt((1/length(sds))*sum(sds^2))
		}
	}
	
	# summarise sd_eta similarly to sd_epsilon
	toRemove2 <- which(is.na(sd_etas))
	if(length(toRemove2)>0){
		sd_eta <- sqrt((1/(l-r))*sum((sd_etas[-toRemove2]^2)*(ljs[-toRemove2]-1)))
	}
	if(length(toRemove2)==0){
		sd_eta <- sqrt((1/(l-r))*sum((sd_etas^2)*(ljs-1)))
	}

	cat("calculated Rocke-Durbin error model parameters \n")
	list(alpha=alpha,sd_epsilon=sd_epsilon,sd_eta=sd_eta)
}


