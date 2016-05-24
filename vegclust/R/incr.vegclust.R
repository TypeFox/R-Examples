incr.vegclust<-function(x, method="NC", ini.fixed.centers = NULL, min.size = 10, alpha = 0.5, nstart=100, fix.previous = TRUE, ...) {
  METHODS <- c("NC", "HNC", "NCdd", "HNCdd")
  method <- match.arg(method, METHODS)  
  k=1
  cat(paste("Initial vegclust with one new group..."))
  vc<-vegclust(x,mobileCenters=k, fixedCenters=ini.fixed.centers, method=method, nstart=nstart,...)
  #Number of seed objects 
  noise<-vc$memb[,ncol(vc$memb)]>alpha
  cat(paste("Number of remaining cluster seeds: ", sum(noise),"\n",sep=""))
  #Stop before continuing if there are not enough seeds object or the first cluster is too small
  if(sum(colSums(vc$memb[,-ncol(vc$memb), drop=FALSE]>alpha)<min.size)>0) {
    cat("The initial cluster was too small.\n")
    return(NULL)
  }
  if(sum(noise)==0) {
    cat("Not enough objects to act as cluster seeds. Stopping.\n")
    return(vc)
  } 
  #Continue
  cont = TRUE
	while(cont) {
	  k = k+1
	  vcold = vc
	  vcbest = NULL
    if(nstart<sum(noise)) seeds = sample(which(noise),nstart)
    else seeds = which(noise)
	  cat(paste("Vegclust with", k,"new groups..."))
    fixed = vcold$mobileCenters[1,]
	  if(!is.null(vcold$fixedCenters)) fixed = rbind(fixed,vcold$fixedCenters)
	  for(i in 1:length(seeds)) {
      if(fix.previous) vc<-vegclust(x,mobileCenters=x[seeds[i],],
                                    fixedCenters=fixed, method=method, nstart=nstart,...)
      else vc<-vegclust(x,mobileCenters=rbind(fixed,x[seeds[i],]), 
                        fixedCenters = ini.fixed.centers, method=method, nstart=nstart,...)
      if(is.null(vcbest)) vcbest = vc
      else if(vc$functional<vcbest$functional) vcbest = vc
    }
    vc = vcbest
	  noise<-vc$memb[,ncol(vc$memb)]>alpha
	  cat(paste("Number of remaining cluster seeds: ", sum(noise),"\n",sep=""))
	  cont = (sum(colSums(vc$memb[,-ncol(vc$memb), drop=FALSE]>alpha)<min.size)==0 && sum(noise)>0)
	}
	if(sum(colSums(vc$memb[,-ncol(vc$memb), drop=FALSE]>alpha)<min.size)>0) {
	  cat("Some of the current clusters are too small. ")
	  cat(paste("Returning vegclust with", k-1,"new group(s).\n"))
	  return(vcold)
	} else if(sum(noise)==0) {
	  cat("Not enough objects to act as cluster seeds. Stopping.\n")
	  return(vc)
	}
}