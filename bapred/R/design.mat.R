design.mat <-
function(mod, numCov){

	tmp <- which(colnames(mod) == 'Batch')
	tmp1 <- as.factor(mod[,tmp])
	# cat("Found",nlevels(tmp1),'batches\n')
	design <- build.design(tmp1,start=1)

	if(!is.null(numCov)) {
		theNumCov = as.matrix(mod[,numCov])
		mod0 = as.matrix(mod[,-c(numCov,tmp)])
	} else 	mod0 = as.matrix(mod[,-tmp])

	ncov <- ncol(mod0)
	
	# cat("Found",ncov,' categorical covariate(s)\n')
	if(!is.null(numCov)) cat("Found",ncol(theNumCov),' continuous covariate(s)\n')
	if(ncov>0){
		for (j in 1:ncov){
			tmp1 <- as.factor(as.matrix(mod0)[,j])
			design <- build.design(tmp1,des=design)
			}
		}
	if(!is.null(numCov)) design = cbind(design,theNumCov)
	return(design)

}
