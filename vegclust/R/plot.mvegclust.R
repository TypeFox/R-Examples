plot.mvegclust<-function(x, type="hnc", excludeFixed=TRUE, verbose=FALSE, ylim=NULL, xlab=NULL, ylab=NULL, maxvar=0.6, minsize=20,...) {
  minclasses=999999
  maxclasses=0
  nplots = NA
  getNclasses<-function(a,excludeFixed=FALSE) {
  	  n = names(a$dist2clusters)
  	  if(excludeFixed) n = n[substr(n,1,1)!="F"]
  	  return(length(n))
  	}
	for(i in 1:length(x)) {
		if(!is.null(x[[i]])) {
			nclasses = getNclasses(x[[i]], excludeFixed)
			maxclasses = max(maxclasses,nclasses)
			minclasses = min(minclasses,nclasses)
			nplots = nrow(x[[i]]$memb)
		}
	}
	if(verbose) print(minclasses)
	if(verbose) print(maxclasses)
	if(verbose) print(nplots)
  if(type=="hnc") {
		nnoise<-rep(NA,maxclasses-minclasses+1)
		for(i in 1:length(x)){
			if(!is.null(x[[i]])) {
				nclasses = getNclasses(x[[i]] , excludeFixed)
				nnoise[nclasses-minclasses+1]<-sum(x[[i]]$memb[,length(names(x[[i]]$memb))]>0.5)
			}
		}
		noisec<-rep(NA,length(nnoise))
		noisec[1] = nnoise[1]-nplots
		for(i in 2:length(nnoise)) {
				noisec[i]<-nnoise[i]-nnoise[i-1]
		}
		plot(minclasses:maxclasses,noisec,type="l",ylab="Changes in the Noise class", xlab="# groups",...)
		axis(1, xaxp=c(minclasses,maxclasses, (maxclasses-minclasses)), tcl = -0.2, labels=FALSE)
		abline(h=0, col="red", lty=3) 
		invisible(noisec)
  	}  else if(type=="hmemb") {
		nnoise<-rep(NA,maxclasses-minclasses+1)
		nass<-rep(NA,maxclasses-minclasses+1)
		nint<-rep(NA,maxclasses-minclasses+1)
		for(i in 1:length(x)){
			if(!is.null(x[[i]])) {
				nclasses = getNclasses(x[[i]] , excludeFixed)
				nass[nclasses-minclasses+1]<-sum(x[[i]]$memb[,1:nclasses]>0.5)
				nnoise[nclasses-minclasses+1]<-sum(x[[i]]$memb[,length(names(x[[i]]$memb))]>0.5)
				nint[nclasses-minclasses+1] = nplots - nass[nclasses-minclasses+1] - nnoise[nclasses-minclasses+1]
			}
		}
		propNoise<-nnoise/nplots
		propInt<-nint/nplots
		propAss<-nass/nplots
		if(is.null(ylim)) ylim <- c(0,1)
		if(is.null(xlab)) xlab <- "Number of groups"
		if(is.null(ylab)) ylab <- "Proportion of all objects"
		plot(minclasses:maxclasses, propNoise,type="n", ylim=ylim,axes=TRUE, frame.plot=FALSE, ylab=ylab, xlab=xlab,...)
		axis(1, xaxp=c(minclasses,maxclasses, (maxclasses-minclasses)), tcl = -0.2, labels=FALSE)
		lines(minclasses:maxclasses, propNoise, lty=1)
		lines(minclasses:maxclasses, propAss, lty=2)
		lines(minclasses:maxclasses, propInt, lty=3)
		legend("topright", legend=c("Noise class","Assigned","Transitional"), lty=c(1,2,3), bty="n")
		invisible(data.frame(propAss,propInt,propNoise))
  	} else if(type=="var") {
		minvar<-rep(NA,maxclasses-minclasses+1)
		meanvar<-rep(NA,maxclasses-minclasses+1)
		maxvar<-rep(NA,maxclasses-minclasses+1)
		for(i in 1:length(x)){
			if(!is.null(x[[i]])) {
				nclasses = getNclasses(x[[i]] , excludeFixed)
				cv<-clustvar(x[[i]])[1:nclasses]
				minvar[nclasses-minclasses+1]<-min(cv)
				meanvar[nclasses-minclasses+1]<-mean(cv)
				maxvar[nclasses-minclasses+1]<-max(cv)
			}
		}
		plot(minclasses:maxclasses,meanvar, ylim=c(min(minvar, na.rm=TRUE),max(maxvar, na.rm=TRUE)),type="l", ylab="Fuzzy cluster variance", xlab="# groups",...)
		axis(1, xaxp=c(minclasses,maxclasses, (maxclasses-minclasses)), tcl = -0.2, labels=FALSE)
		lines(minclasses:maxclasses, minvar, lty=2)
		lines(minclasses:maxclasses, maxvar, lty=2)
		invisible(data.frame(minvar,meanvar,maxvar))
  	} else if(type=="hcs") {
		mincs<-rep(NA,maxclasses-minclasses+1)
		meancs<-rep(NA,maxclasses-minclasses+1)
		maxcs<-rep(NA,maxclasses-minclasses+1)
		for(i in 1:length(x)){
			if(!is.null(x[[i]])) {
				nclasses = getNclasses(x[[i]] , excludeFixed)
				cs<-colSums(x[[i]]$memb>0.5)[1:nclasses]
				mincs[nclasses-minclasses+1]<-min(cs)
				meancs[nclasses-minclasses+1]<-mean(cs)
				maxcs[nclasses-minclasses+1]<-max(cs)
			}
		}
		if(is.null(ylim)) ylim <- c(min(mincs, na.rm=TRUE),max(maxcs, na.rm=TRUE))
		if(is.null(xlab)) xlab <- "Number of groups"
		if(is.null(ylab)) ylab <- "(Hard) cluster size"
		plot(minclasses:maxclasses,meancs, ylim=ylim,type="l", ylab=ylab, xlab=xlab,...)
		axis(1, xaxp=c(minclasses,maxclasses, (maxclasses-minclasses)), tcl = -0.2, labels=FALSE)
		lines(minclasses:maxclasses, mincs, lty=2)
		lines(minclasses:maxclasses, maxcs, lty=2)
		invisible(data.frame(mincs,meancs,maxcs))
  	} else if(type=="valid") {
		lowcs<-rep(NA,maxclasses-minclasses+1)
		highvar<-rep(NA,maxclasses-minclasses+1)
		invalid<-rep(NA,maxclasses-minclasses+1)
		for(i in 1:length(x)){
			if(!is.null(x[[i]])) {
				nclasses = getNclasses(x[[i]] , excludeFixed)
				cs<-colSums(x[[i]]$memb>0.5)[1:nclasses]
				cv<-clustvar(x[[i]])[1:nclasses]
				lowcs[nclasses-minclasses+1]<-sum(cs<minsize)
				highvar[nclasses-minclasses+1]<-sum(cv> maxvar)
				invalid[nclasses-minclasses+1]<-sum(cv> maxvar & cs<minsize)
			}
		}
		if(is.null(ylim)) ylim <- c(0,max(c(invalid, highvar, lowcs)+1,na.rm=TRUE))
		if(is.null(xlab)) xlab <- "Number of groups"
		if(is.null(ylab)) ylab <- "Number of invalid groups"
		plot(minclasses:maxclasses,invalid, ylim=ylim,type="l", ylab=ylab, xlab=xlab,...)
		axis(1, xaxp=c(minclasses,maxclasses, (maxclasses-minclasses)), tcl = -0.2, labels=FALSE)
		axis(2, yaxp=c(ylim[1],ylim[2], ylim[2]-ylim[1]), tcl = -0.2, labels=FALSE)
		lines(minclasses:maxclasses, highvar, lty=2)
		lines(minclasses:maxclasses, lowcs, lty=3)
		legend("topleft",legend=c("High variance","Low membership", "High variance & low membership"), lty=c(2,3,1), bty="n")
		invisible(data.frame(lowcs,highvar,invalid))
  	}
}