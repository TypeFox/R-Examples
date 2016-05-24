tunelearnPattern <- function(x, y, unlabeledx=NULL, nfolds=5, segmentlevels=c(0.25,0.5,0.75), random.split=0, 
						mindepth=4, maxdepth=8, depthstep=2, ntreeTry=25, target.diff=TRUE, segment.diff=TRUE, ...) {

	ndepthlev <- 1+(maxdepth-mindepth)/depthstep
	
	error_rates <- matrix(0,length(segmentlevels)*ndepthlev,nfolds)
	param_combination <- matrix(0,length(segmentlevels)*ndepthlev,2)

	noftrain_labeled <- nrow(x)

	#cross-validation indices
	id <- sample(rep(seq.int(nfolds), length.out=noftrain_labeled))
    cv <- lapply(seq.int(nfolds), function(x) list(
         train=which(x!=id),
         test=which(x==id)
    ))
    
	for(fold in 1:nfolds){
		trainind <- cv[[fold]]$train
		testind <- cv[[fold]]$test
		ntrain <- length(trainind); classtr=y[trainind];
		ntest <- length(testind); classtst=y[testind];
		
		cnt <- 1
		for(l in 1:length(segmentlevels)){
			if(!is.null(unlabeledx)){
				ensemble <- learnPattern(rbind(x[trainind,],unlabeledx),segment.factor=segmentlevels[l],random.seg=random.split,
					target.diff=target.diff,segment.diff=segment.diff,ntree=ntreeTry,maxdepth=maxdepth)	
			} else {
				ensemble <- learnPattern(x[trainind,],segment.factor=segmentlevels[l],random.split=random.split,
					target.diff=target.diff,segment.diff=segment.diff,ntree=ntreeTry,maxdepth=maxdepth)		
			}
			
			tempdepth <- mindepth
			while(tempdepth<=maxdepth){	
				sim <- computeSimilarity(ensemble,x[testind,],x[trainind,],maxdepth=tempdepth)		
				id <- apply(sim,1,which.min)					
				predicted <- classtr[id]
				error_rates[cnt,fold] <- 1-sum(classtst==predicted)/ntest	
				cnt <- cnt+1
				tempdepth <- tempdepth+depthstep
			}
			rm(sim,ensemble); gc();
		}
	} 
	
	#parameter combinations
	cnt <- 1
	for(l in 1:length(segmentlevels)){
		tempdepth <- mindepth
		while(tempdepth<=maxdepth){
			param_combination[cnt,1:2] <- c(segmentlevels[l],tempdepth)	
			tempdepth <- tempdepth+depthstep
			cnt <- cnt+1
		}
	}
	
	avg_error <- apply(error_rates,1,mean)
	minindex <- which(avg_error==min(avg_error))
	best_depth <- param_combination[minindex[1],2]
	best_seg <- param_combination[minindex[1],1]
			
	res <- list(params=param_combination,errors=error_rates,
				best.error=min(avg_error),best.seg=best_seg,best.depth=best_depth,random.split=random.split)   	
}
