GESTr <- function(x,dist.theta=0.05,merge.overlap=0.1,verbose=FALSE){

	RDparameters <- fitRockeDurbin(x,theta=dist.theta)

	epsilons <- rnorm(1000,mean=0,sd=RDparameters$sd_epsilon)
        epsilons <- (epsilons-mean(epsilons))/sd(epsilons)
        etas <- rnorm(1000,mean=0,sd=RDparameters$sd_eta)
        etas <- (etas-mean(etas))/sd(etas)

	confidences_lo <- array(0,dim=dim(x))
	confidences_hi <- array(0,dim=dim(x))
	
	for(i in c(1:nrow(x))){
		vals <- x[i,]
		
		confidence_lo <- rep(0,length(vals))
	        confidence_hi <- rep(0,length(vals))

		geneModel <- fitGMM(vals,RDparameters)
		geneModel$G <- max(geneModel$classification)
	        mergeList <- mergeComponents(geneModel,overlap=merge.overlap)

		if(geneModel$G > 2){
                	allClasses <- 1:geneModel$G
	                intermediateClasses <- setdiff(allClasses,c(mergeList[[1]],mergeList[[2]]))
	                if(length(intermediateClasses) > 0){
	                        intermediatevals <- list()
	                        for(k in 1:length(intermediateClasses)){
	                                intermediatevals[[k]] <- rep(0,length(vals))
                	        }
        	        }
	        }

		if(geneModel$G>0){
	                if(verbose) cat(paste("calculating confidence values for each sample, for gene",i,"with",geneModel$G,"model components \n"))
	                for(j in 1:length(vals)){
	                        y <- 2^vals[j]
	                        if(y < RDparameters$alpha){
	                                # if value is less than background, just use noise
	                                simulatedVals <- log(epsilons/exp(etas),base=2)
	                                simulatedVals <- simulatedVals[!is.na(simulatedVals)]
	                        }
	                        else{
	                                simulatedVals <- log((y-(RDparameters$alpha + epsilons))/exp(etas),base=2)
	                                simulatedVals <- simulatedVals[!is.na(simulatedVals)]
	                        }
	                        if(geneModel$G==1){
	                                confidence_lo[j] <- mean(unlist(lapply(simulatedVals,function(x,mu,sd){1-pnorm((x-mu)/sd)},mu=geneModel$parameters$mean,sd=sqrt(geneModel$parameters$variance$sigmasq))))
        	                        confidence_hi[j] <- mean(unlist(lapply(simulatedVals,function(x,mu,sd){pnorm((x-mu)/sd)},mu=geneModel$parameters$mean,sd=sqrt(geneModel$parameters$variance$sigmasq))))
	                        }
	                        if(geneModel$G==2){
	                                zs <- estep("V",data=simulatedVals,parameters=geneModel$parameters)
	                                confidence_lo[j] <- mean(zs$z[,1])
	                                confidence_hi[j] <- mean(zs$z[,2])
	                        }
	                        if(geneModel$G > 2){
	                                zs <- estep("V",data=simulatedVals,parameters=geneModel$parameters)
	                                if(length(mergeList[[1]])>1){
	                                        confidence_lo[j] <- length(mergeList[[1]])*mean(zs$z[,mergeList[[1]]])
	                                }
	                                else{
	                                        confidence_lo[j] <- mean(zs$z[,1])
	                                }
	                                if(length(mergeList[[2]])>1){
	                                        confidence_hi[j] <- length(mergeList[[2]])*mean(zs$z[,mergeList[[2]]])
	                                }
	                                else{
	                                        confidence_hi[j] <- mean(zs$z[,geneModel$G])
	                                }
	                                if(length(intermediateClasses) > 0){
	                                        for(k in 1:length(intermediateClasses)){
	                                                intermediatevals[[k]][j] <- mean(zs$z[,intermediateClasses[k]])
	                                        }
	                                }
	                        }
	                }
		}

		if(verbose) cat("ensuring monotonicity of confidence-assignments \n")
	        if(geneModel$G < 3){
	                classScoresList <- list(confidence_lo,confidence_hi)
	        }
	        else{
	                classScoresList <- list()
	                classScoresList[[1]] <- confidence_lo
	                if(length(intermediateClasses) > 0){
	                        for(k in 1:length(intermediateClasses)){
	                                classScoresList[[length(classScoresList)+1]] <- intermediatevals[[k]]
	                        }
	                }
	                classScoresList[[length(classScoresList)+1]] <- confidence_hi
	        }
	
        	# sort confidences on expression order
	        sortOrder <- order(vals,decreasing=FALSE)
	        for(k in 1:length(classScoresList)){
	                classScoresList[[k]] <- classScoresList[[k]][sortOrder]
	        }
	
	        # ensure monotonicity of hi and lo confidence scores
	        if(geneModel$G < 3){
	                monotonicConfidences <- getMonotonicConfidences(classScoresList)
	        }
	        else{
	                monotonicConfidences <- getMonotonicConfidences_multiclass(classScoresList,vals[sortOrder],verbose)
	        }
	        confidence_lo[sortOrder] <- monotonicConfidences$confidences_lo
	        confidence_hi[sortOrder] <- monotonicConfidences$confidences_hi
	
		confidences_lo[i,] <- confidence_lo
		confidences_hi[i,] <- confidence_hi
	}
	tx.matrix <- (1+confidences_hi-confidences_lo)/2
	rownames(tx.matrix) <- rownames(x)
	colnames(tx.matrix) <- colnames(x)
	tx.matrix
}
