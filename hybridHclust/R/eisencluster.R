eisenCluster <- function(x,method,compatible=TRUE,verbose=FALSE){

cordist <- function(x,absolute=FALSE) {
	cd <- cor(t(x),use='pairwise.complete.obs')
	if(absolute) return(1-abs(cd)) else return(1-cd)
}

uncent.cor2dist <- function(x,i,absolute=FALSE){
# Calculates uncentered correlation distance
# between row i of x and all other rows of x.
	myfun <- function(y){
		x1 <- x[i,]
		x1[is.na(y)] <- NA
		y[is.na(x1)] <- NA
		sum(y*x1/sqrt(sum(y*y,na.rm=TRUE)*sum(x1*x1,na.rm=TRUE)),na.rm=TRUE)
		#mean(y*x1/sqrt(mean(y*y,na.rm=TRUE)*mean(x1*x1,na.rm=TRUE)),na.rm=TRUE)
	}
	if(absolute) return(1-abs(apply(x,1,myfun))) else return(1-apply(x,1,myfun))
}


uncent.cordist <- function(x,absolute=FALSE) {
# calculates all pairwise distances using correlation distance
	cd <- matrix(0,nrow(x),nrow(x))
	for (i in 1:nrow(x)) {
		cd[,i] <- uncent.cor2dist(x,i,absolute)
	}
	#browser()
	return(cd)  # value of "absolute" already used above
}



#---------------------------------------------------------------------- 

	METHODS <- c('euclidean','squared.euclidean','correlation','uncentered.correlation')
	method <- pmatch(method,METHODS)
	if (method==1) {
		dfun <- function(x,absolute=FALSE){as.matrix(daisy(x,metric='euclidean'))} 
		d2fun <- function(x,i,absolute=FALSE) {
			diff <- (x - x[rep(i,nrow(x)),])^2
			multiplier <- ncol(x)/apply(!is.na(diff),1,sum)
			sqrt(apply(diff,1,sum,na.rm=TRUE)*multiplier)
		}
	} else if (method==2) {
		dfun <- function(x,absolute=FALSE){as.matrix(daisy(x,metric='euclidean')^2)} 
		d2fun <- function(x,i,absolute=FALSE) {
			diff <- (x - x[rep(i,nrow(x)),])^2
			multiplier <- ncol(x)/apply(!is.na(diff),1,sum)
			apply(diff,1,sum,na.rm=TRUE)*multiplier
		}
	} else if (method==3) {
		dfun <- function(x,absolute=FALSE) {cordist(x,absolute)}
		d2fun <- function(x,i,absolute=FALSE){
		# need to be careful about missing values and scaling.  Corresponds
		# to use='pairwise.complete.obs'
			x <- x[,!is.na(x[i,])]
			x1 <- t(apply(x,1,function(x) {(x-mean(x,na.rm=TRUE))/sqrt(var(x,na.rm=TRUE))}))
			x2 <- x[rep(i,nrow(x)),]
			x2[is.na(x)] <- NA
			x2 <- t(apply(x2,1,
				function(x) {(x-mean(x,na.rm=TRUE))/sqrt(var(x,na.rm=TRUE))}))
			nob <- apply(!is.na(x),1,sum)
			cd <- apply(x1*x2,1,sum,na.rm=TRUE)/(nob-1)
			if (absolute) return(1-abs(cd)) else return(1-cd)
		}
	} else if (method==4) {
		d2fun <- function(x,i,absolute=FALSE){uncent.cor2dist(x,i,absolute)}
		dfun <- function(x,absolute=FALSE) {uncent.cordist(x,absolute)}
	} else stop(paste('Method must be one of',METHODS,'\n'))

	if ((nrow(x)==ncol(x)) && (all(x==t(x))))
		cat('Warning: x looks like it might be a distance matrix\n')

	descendants <- function(m,k){
	# the done object indicates what rows of m were used
		done <- k
		if (m[k,1]<0) left <- -m[k,1]
		else {
			junk <- descendants(m,m[k,1])
			left <- junk[[1]]
			done <- c(done,junk[[2]])
		}
		if (m[k,2]<0) right <- -m[k,2]
		else {
			junk <- descendants(m,m[k,2])
			right <- junk[[1]]
			done <- c(done,junk[[2]])
		}
		return(list(c(left,right),done))
	}

	findminsym <- function(m,include.diag=FALSE){
	# Returns the indices of elements of m that are minimum.
	# Each row of matrix returned gives (row, column)
	# Note that m is symmetric
		if (!include.diag) diag(m) <- 1e10
		w <- which(m==min(m))
		ww <- matrix(0,length(w),2)
		for (i in 1:length(w)) 
			ww[i,] <- c((w[i]-1)%/%nrow(m)+1,w[i]%%nrow(m))
		ww[ww[,2]==0,2] <- nrow(m)
		for (i in 1:nrow(ww)) ww[i,] <- sort(ww[i,])
		return(ww[!duplicated(apply(ww,1,paste,collapse=' ')),, drop=FALSE])
	}

	findminsym <- function(m,include.diag=FALSE){
		if (!include.diag) diag(m) <- 1e10
		m[upper.tri(m)] <- 1e10
		which(m==min(m),arr.ind=TRUE)
	}

	height <- rep(0,nrow(x)-1)
	merge <- matrix(0,nrow(x)-1,2)
	group <- vector('list',nrow(x))
	for (i in 1:nrow(x)) group[[i]] <- list(id=i,lastmerge=0)
	dmat <- dfun(x)
	diag(dmat) <- 1e10
	dmat[is.na(dmat)] <- 1e10

	xorig <- x

	for (k in 1:length(height))
	{
		#print(str(group))
		#print(dmat)
		#print(x)
		#browser()
		if (verbose) cat(k)
		mi <- findminsym(dmat)  # mi short for minindex
		if (nrow(mi)>1) {
			mi <- mi[1,,drop=FALSE]
			cat('\nTied minimum distance occured.  Taking first occurrance\n')
		}
		mi <- sort(mi[1,])
		j1 <- ifelse(group[[mi[1]]]$lastmerge==0,-group[[mi[1]]]$id,
			group[[mi[1]]]$lastmerge)
		j2 <- ifelse(group[[mi[2]]]$lastmerge==0,-group[[mi[2]]]$id,
			group[[mi[2]]]$lastmerge)
		merge[k,] <- c(j1,j2)
		height[k] <- min(dmat)
		keep <- (1:length(group))[-mi[2]]

		if (!compatible){
			# next line is the `right' way to do this, but doesn't give
			# Eisen's results.
			x[mi[1],] <- apply(xorig[c(group[[mi[1]]]$id,group[[mi[2]]]$id),],2,
				mean,na.rm=TRUE)
		} else {
			# for compatability with Eisen's code
			n1 <- length(group[[mi[1]]]$id)
			n2 <- length(group[[mi[2]]]$id)
			x[mi[1],] <- apply(x[rep(mi,c(n1,n2)),],2,mean,na.rm=TRUE)
			#cat('new x has mean',mean(x[mi[1],],na.rm=TRUE),'and sd',
			#	sd(x[mi[1],],na.rm=TRUE),'\n')
			#browser()
		}
		group[[mi[1]]]$id <- c(group[[ mi[1] ]]$id,group[[ mi[2] ]]$id)
		group[[mi[1]]]$lastmerge <- k

		dmat <- dmat[keep,keep,drop=FALSE]
		x <- x[keep,,drop=FALSE]
		#browser()
		if (k<length(height)) {
			# since this is calculating distances between remaining points after
			# merge, we don't need to update it in last loop.
			dmat[mi[1],] <- dmat[,mi[1]] <- d2fun(x,mi[1])
			dmat[mi[1],mi[1]] <- 1e10
		}
		dmat[is.na(dmat)] <- 1e10

		group[[mi[2]]] <- NULL
	}

	storage.mode(merge) <- 'integer'
	myorder <- descendants(merge,length(height))[[1]]
	storage.mode(myorder) <- 'integer'
	realheight <- height
	if (any(diff(height)<0)){
		cat('\nNon-monotone dendogram.  Forcing heights to be monotone\n')
		for (i in 2:length(height))
			if (height[i]<height[i-1]) height[i] <- height[i-1]
	}
	tree <- list(merge=merge,height=height,realheight=realheight,
		order=myorder,labels=NULL,method='average',
		call=match.call())
	class(tree) <- "hclust"
	tree
}

