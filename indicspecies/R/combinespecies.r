combinespecies<-function(X, max.order = 3, min.occ = 0, FUN = min, verbose=FALSE, add.names=TRUE, ...) {
	nsites = nrow(X)
  	spplist = names(X)
  	nspecies = length(spplist)
  	FUN = match.fun(FUN)  	
  	if(verbose) cat(paste("Number of candidate species: ",nspecies,"\n", sep=""))
  	if(length(spplist)==1) stop("At least two species are necessary.")
	if(verbose) cat(paste("Number of sites:",nsites,"\n"))
	order  = min(max.order,nspecies)
	if(verbose) cat(paste("Maximum order of combinations:",order,"\n"))
	totco = 0
	for(j in 1:order) totco = totco + choose(nspecies,j)
	if(verbose) cat(paste("Total number of combinations:",totco,"\n"))
    C = matrix(0, nrow=totco, ncol=nspecies)
    combi = 1
    for(j in 1:min(max.order, nspecies)) {
      co <- combn(nspecies,j)
      for(i in 1:ncol(co)) {
      	C[combi,co[,i]] = 1
      	combi = combi+1
      }
	}
    XC = matrix(0, nrow=nsites, ncol=totco)
    for(r in 1:totco) {
  		if(sum(C[r,])==1) XC[,r]<-X[,C[r,]==1]
    	else XC[,r]<-apply(X[,C[r,]==1],1,FUN=FUN,...)
    }
    if(min.occ>0) {
    	toRem = colSums(ifelse(XC>0,1,0))<min.occ
    	XC = XC[,!toRem]
    	C = C[!toRem,]
		if(verbose) cat(paste("Number of combinations with low occurrence:",sum(toRem),"\n"))
		if(verbose) cat(paste("Final number of combinations:",sum(!toRem),"\n"))
    }
    XC = as.data.frame(XC)
    row.names(XC) = row.names(X)
    if(add.names) {
    	for(i in 1:ncol(XC)) names(XC)[i] = paste(spplist[C[i,]==1], collapse="+")
    } else names(XC) = 1:ncol(XC)
    C = as.data.frame(C)
    names(C) = spplist
    row.names(C) = names(XC)
    return(list(XC = XC, C=C))
}