relate.levels<-function(lower, upper, defuzzify=FALSE, excludeFixed = FALSE, verbose=FALSE, ...) {
  minupperclasses =999999
  maxupperclasses =0
	minlowerclasses=999999
  	maxlowerclasses =0
	for(i in 1:length(upper)) {
		if(!is.null(upper[[i]])) {
			nclasses = length(names(upper[[i]]$dist2clusters))
			maxupperclasses = max(maxupperclasses,nclasses)
			minupperclasses = min(minupperclasses,nclasses)
		}
	}
	for(i in 1:length(lower)) {
		if(!is.null(lower[[i]])) {
			nclasses = length(names(lower[[i]]$dist2clusters))
			minlowerclasses = min(minlowerclasses,nclasses)
			maxlowerclasses = max(maxlowerclasses,nclasses)
		}
	}
	#Prepare output matrices
	maxnoise = data.frame(matrix(NA,length(minupperclasses:maxupperclasses), length(minlowerclasses:maxlowerclasses)))
	row.names(maxnoise) = minupperclasses:maxupperclasses
	names(maxnoise) = minlowerclasses:maxlowerclasses
	minmaxall = maxnoise
	nnoise = maxnoise
	minallsize = maxnoise
	empty = maxnoise
	
	#Loop over lower clustering
	for(j in 1:length(lower)) {
		if(!is.null(lower[[j]])) {
			nlowerclasses = length(names(lower[[j]]$dist2clusters))
			colind<-nlowerclasses-minlowerclasses+1
			if(excludeFixed) {
			  nfixed = sum(substr(names(lower[[j]]$dist2clusters),1,1)=="F")
        nlowerclasses = nlowerclasses - nfixed
        colind<-nlowerclasses-minlowerclasses+1 + nfixed
			}
			if(verbose) cat(paste("Number of lower classes:", nlowerclasses,"\n"))
			#Loop over upper clustering
			for(i in 1:length(upper)) {
				if(!is.null(upper[[i]])){
					if(verbose) cat(".")
					nupperclasses = length(names(upper[[i]]$dist2clusters))
					rowind<-nupperclasses-minupperclasses+1
					if(excludeFixed) {
            nfixed = sum(substr(names(upper[[i]]$dist2clusters),1,1)=="F")
            nupperclasses = nupperclasses - nfixed
            rowind<-nupperclasses-minupperclasses+1 + nfixed
					}
					memb<-crossmemb(lower[[j]], upper[[i]])
					if(lower[[j]]$method=="NC") memb = memb[-nrow(memb),]
					if(defuzzify) memb <- defuzzify(memb, ...)$memb
					minmaxall[rowind,colind] = min(apply(memb[1:nlowerclasses,1:nupperclasses],2,max))
					minallsize[rowind,colind] = min(apply(memb[1:nlowerclasses,1:nupperclasses],2,sum))
					empty[rowind,colind] = sum(colSums(defuzzify(memb[1:nlowerclasses,1:nupperclasses], method="max")$memb)==0)
					#print(empty[rowind,colind])
					if(upper[[i]]$method=="NC") {
						maxnoise[rowind,colind]=max(memb[1:nlowerclasses,ncol(memb)])
						nplots = colSums(defuzzify(lower[[j]], method="cut", alpha=0.5)$memb)
						nnoise[rowind,colind]=sum(memb[1:nlowerclasses,ncol(memb)]>0.5)
					 	#print(nnoise[rowind,colind])
					}
				}
			}
		}
		if(verbose) cat("done.\n")
	}	

	return(list(nnoise =nnoise, maxnoise=maxnoise, minallsize=minallsize, minmaxall = minmaxall, empty= empty))
}