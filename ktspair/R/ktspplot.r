ktspplot <- function(ktspobj, select=NULL){
	grp <- ktspobj$grp
	labels <- ktspobj$labels
	dat <- ktspobj$ktspdat
	index <- ktspobj$index
	k <- ktspobj$k

	if(is.null(select)){
		cat("Number of TSPs: ",k,"\n")
		for(i in 1:k){
			par(mar=c(4,4,4,4))
			plot(dat[i,],dat[(i + k),],xlab=paste("Gene:",rownames(dat)[i],"Expression"),ylab=paste("Gene:",rownames(dat)[(i+k)],"Expression"),type="n")
			mtext(paste("Groups:",labels[1],"= Red |",labels[2],"= Blue; Score:",round(ktspobj$ktspscore[i],3)),line=1)
			points(dat[i,grp==0],dat[(i+k),grp==0],col="red",pch=19)
			points(dat[i,grp==1],dat[(i+k),grp==1],col="blue",pch=19)
			abline(c(0,1),lwd=2)
			readline(paste("TSP",i,": Hit return for next TSP.\n"))
		}
	}

	if(!is.null(select) && (select <1 || select > k)){stop("The selected pair of genes is not available")}

	if(!is.null(select)){
		par(mar=c(4,4,4,4))
		plot(dat[select,],dat[(select + k),],xlab=paste("Gene:",rownames(dat)[select],"Expression"),ylab=paste("Gene:",rownames(dat)[(select+k)],"Expression"),type="n")
		mtext(paste("Groups:",labels[1],"= Red |",labels[2],"= Blue; Score:",round(ktspobj$ktspscore[select],3)),line=1)
		points(dat[select,grp==0],dat[(select+k),grp==0],col="red",pch=19)
		points(dat[select,grp==1],dat[(select+k),grp==1],col="blue",pch=19)
		abline(c(0,1),lwd=2)   
	}
}

