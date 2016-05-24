PoisDev<-function(counts,design,log.offset,print.progress=TRUE){
n<-ncol(counts)
### For one factor designs, MLE's for Poisson GLM have simple closed form, so we can avoid one-gene-at-a-time GLM fitting.
	if(is.vector(design)) {

		##offset should NOT be given on log scale here
		offset<-rep(1,length(design)); if(!is.null(log.offset)) offset<-exp(log.offset)

		counts<-as.matrix(counts)
		means<-counts

		parms<-NULL
		for(i in 1:length(unique(design))){
		if(sum( design==unique(design)[i])>1) means[,design==unique(design)[i]]<-rowSums(counts[,design==unique(design)[i]])/sum(offset[design==unique(design)[i]])
		if(sum( design==unique(design)[i])==1) means[,design==unique(design)[i]]<-counts[,design==unique(design)[i]]/offset[design==unique(design)[i]]
		parms<-cbind(parms,means[,design==unique(design)[i]][,1])
		}
		parms<-cbind(log(parms[,1]),log(parms[,-1]/parms[,1]))
		colnames(parms)<-c("(Intercept)",unique(design)[-1])
		means<-t(t(means)*offset)

		### 0's require special attention since 0^0=1, but 0*log(0)=NaN
		deviance<-means-counts
		deviance[counts!=0]<-deviance[counts!=0]+(counts*log(counts/means))[counts!=0]
		deviance<-2*rowSums(deviance)
	}

	if(is.matrix(design)) {
	### For multi-factor designs, the first column of each element (matrix) of design.list should be a column of 1's, pertaining to the intercept. 

		deviance<-rep(NA,nrow(counts)); means<-matrix(NA,nrow(counts),ncol(counts)); parms<-matrix(NA,nrow(counts),ncol(design))

		### For each gene and given design matrix, fit GLM to find model parameters (for mean structure) that optimize quasi-likelihood
		for(i in 1:nrow(counts)){
			### If wanted, provide running progress update (eventually once every 5000 genes) 
			if(i%in%c(2,10,100,500,1000,2500,4000,5000*(1:200))&print.progress) print(paste("Analyzing Gene #",i))
			
			### Fit GLM
			res<-withCallingHandlers(
				glm(counts[i,]~design[,-1],family="poisson",offset=log.offset,method=glm.fit3)	##offset should be given on log scale here
				, simpleWarning=ignorableWarnings
			)
			
			### Save optimized means (used in Pearson's dispersion estimator)
			means[i,]<-res$fitted.values
			parms[i,]<-res$coefficients

			### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
			deviance[i]<-res$deviance
		}
	}

return(list(dev=deviance,means=means,parms=parms))
}