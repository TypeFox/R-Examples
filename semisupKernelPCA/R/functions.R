# perform the kernel PCA projection step from an existing kernel matrix
computeProjectionFromKernel <- function(kernel, dims=2, eigentype=c("basic", "irlba")) {
	eigentype <- match.arg(eigentype)
	eigenfunc <- switch(eigentype, basic=function(x) eigen(x, symmetric=TRUE), irlba=function(x) irlba(x, nu=dims, nv=0, adjust=dims))
	nelts <- dim(kernel)[1]
	if (dims > nelts) stop("dims cannot be greater than data rank")
	ones <- matrix(1/nelts, nrow=nelts, ncol=nelts)
	kernel <- kernel - ones%*%kernel - kernel%*%ones + ones %*%kernel%*%ones
	eigendecomp <- eigenfunc(kernel)
	members <- switch(eigentype, basic=c("values", "vectors"), irlba=c("d", "u"))
	# restric processing to dimensions actually used (ie dims)
	values <- eigendecomp[[members[1]]][1:dims]
	vectors <- eigendecomp[[members[2]]][,1:dims,drop=FALSE]
	# normalize eigenvectors according to (12.81)
	vectors <- sweep(vectors, 2, sqrt(values*nelts), "/") 
	# then compute projections according to (12.82), in matrix form
	projection <- kernel %*% vectors # transposed version of Bishop's description
}




computeKernel <- function(data, type=c("gaussian", "pgaussian")) {
	type <- match.arg(type)
	nelts <- dim(data)[1]
	d <- dim(data)[2]
	Ksim <- matrix(0,nrow=nelts, ncol=nelts)
	# first, standardize data
	means <- apply(data[,1:d], 2, mean)
	sds <- apply(data[,1:d], 2, sd)
	sds[sds==0] <- 1 # manage NULL case
	data[,1:d] <- sweep(data[,1:d], 2, means, "-")
	data[,1:d] <- sweep(data[,1:d], 2, sds, "/")

	# compute p and sigma parameters, that will be used later for the p-gaussian kernel	
	# (see Francois2005 for details)
	# NB : D_N = 5% quantile, and D_F = 95% quantile (necessary for p to remain positive)
	# from there, expressions for sigma seem to have been swapped (D_N with 0.95 and reciprocally)
	# => add as a footnote
	Ksim <- dist(data[,1:d]) # contains all unique non-trivial pairwise distances (n.(n-1)/2)
	if (type == "pgaussian") {
		quants <- quantile(Ksim, probs=c(0.05, 0.95))
		p <- log(log(0.05)/log(0.95)) / log(quants[2] / quants[1])
		sigma <- quants[2] / (-log(0.05))^(1/p)	
	} else if (type == "gaussian") {
		p <- 2
		sigma <- max(Ksim)
	} else {
		stop("unsupported kernel type")
	}

	# cast to sim matrix
	Ksim <- as.matrix(Ksim)

	# change distance to similarity values
	if (any(type==c("gaussian", "pgaussian"))) {
		Ksim <- exp((-Ksim^p) / sigma^p)
	}

	return(Ksim)

}


linkinfo <- function(mat) {
	res <- list()
	n <- dim(mat)[1]
	# use binary matrices for representation (should more easily match the kernel representation)
	# makes sense as we use binary relations.
	res$link <- matrix(FALSE, nrow=n, ncol=n)
	res$notlink <- matrix(FALSE, nrow=n, ncol=n)
	class(res) <- "linkinfo"
	return(res)
}

length.linkinfo <- function(x, ...) {
	return(dim(x$link)[1])
}

add <- function(object, ...) {
	UseMethod("add")
}

add.linkinfo <- function(object, inda, indb, type=c("link", "notlink"), ...) {
	type <- match.arg(type)
	n <- length(object)
	# support multiple args	
	
	# is the couple in the correct range ? single numeric values ?
	if (any(inda < 1) || any(inda > n) || any(indb < 1) || any(indb > n)) {
		stop("inda and indb must be in the range of the linkinfo object")
	}
		
	if (length(inda) != length(indb)) {
		stop("inda and indb must have same length")
	}
	
	if (length(inda) == 0) {
		stop("at least one couple must be specified")
	}
	
	if (any(inda == indb)) {
		stop("use with inda != indb")
	}
	
	# shape as matrix for indexing
	inds <- cbind(inda, indb)
	revinds <- cbind(indb, inda)
	
	# does the couple already exist ?
	if (any(object$link[inds]) || any(object$notlink[inds])) {
		stop("[inda, indb] contains at least 1 registered couple")
	}	
	
	# register the couple
	if (type == "link") {
		object$link[inds] <- object$link[revinds] <- TRUE
	} else {
		object$notlink[inds] <- object$notlink[revinds] <- TRUE
	}
	
	return(object)
}

# check consistency ?
# - if 1 is linked to 2, and 2 is linked to 3, it is reasonable that 1 is linked to 3
# - if 1 is linked to 2, and 2 is not linked to 3, it is reasonable that 1 is not linked to 3
# - no other constraint : indeed, if 1 is not linked to 2, and 2 is not linked to 3, nothing can be concluded about 1 and 3 -> they may as well be linked. If they are, contraint 2 is respected.
# may be understood under terms of relations : L is transitive, and we have to respect the logical induction :
# U o L(1,3) -> U(1,3)
check <- function(object, ...) {
	UseMethod("check")		
}


check.linkinfo <- function(object, ...) {
	templink <- transform(object)
	
	# consistency necessary and sufficient condition : matrices/relations cannot "overlap"
	if (any(templink$link & templink$notlink)) {
		return(FALSE) 
	} else {
		return(TRUE)
	}
}

transform.linkinfo <- function(`_data`, ...) {
	# apply transitivity of relation 1
	trans1 <- function(vec) {
		vec <- which(vec==TRUE) # convert to indices
		if (length(vec) > 1) {
			inds <- cbind(rep(vec, times=length(vec)), rep(vec, each=length(vec))) # convert to index matrix
			templink$link[inds] <<- TRUE
		}
	}
	
	# apply transitivity of relation 2
	trans2 <- function(vecN, vecL) {
		if ((length(vecN) > 0) && (length(vecL) > 0)) {
			#vec <- c(vecN, vecL) # no : asymetry must be respected.
			inds <- cbind(rep(vecL, times=length(vecN)), rep(vecN, each=length(vecL)))
			inds <- rbind(inds, cbind(rep(vecN, each=length(vecL)), rep(vecL, times=length(vecN))))
			templink$notlink[inds] <<- TRUE
		}
	}
	
	# for mapply : get list of cols indexes which are at 1 for each row
	indexes <- function(vec) {
		vec <- which(vec==TRUE)
		return(vec)
	}
	
	# copy original linkinfo
	templink <- `_data`
	
	# trans1 : not interested in the result, linkinfo copy is modified as side effect.
	apply(`_data`$link, 1, trans1)
	diag(templink$link) <- FALSE # reset diagonal before trans2 (side effect for simpler algorithm)
	
	# trans2 : first, index lists are extracted. then same remark as for trans1.
	indsL <- lapply(seq(nrow(`_data`$link)), function(i) indexes(`_data`$link[i,]))
	indsN <- lapply(seq(nrow(`_data`$notlink)), function(i) indexes(`_data`$notlink[i,]))
	mapply(trans2, vecN=indsN, vecL=indsL) 
	
	# modification to notlinkset may influence trans2 itself : second application is necessary for stable solution
	#diag(linkinfo$notlink) <- FALSE
	#indsN <- apply(linkinfo$notlink, 1, indexes)
	#mapply(trans2, vecN=indsN, vecL=indsL)
	
	diag(templink$notlink) <- FALSE # reset diagonal
	
	return(templink)
}



del <- function(object, ...) {
	UseMethod("del")		
}


del.linkinfo <- function(object, inda, indb, ...) {
	n <- length(object)
	
	# is the couple in the correct range ? single numeric values ?
	if (any(inda < 1) || any(inda > n) || any(indb < 1) || any(indb > n)) {
		stop("inda and indb must be in the range of the linkinfo object")
	}
	
	if (length(inda) != length(indb)) {
		stop("inda and indb must have same length")
	}

	if (length(inda) == 0) {
		stop("at least one couple must be specified")
	}
	
	if (any(inda == indb)) {
		stop("use with inda != indb")
	}

	# shape as matrix for indexing
	inds <- cbind(inda, indb)
	revinds <- cbind(indb, inda)

	# remove from both linking functions
	object$link[inds] <- object$link[revinds] <- FALSE
	object$notlink[inds] <- object$notlink[revinds] <- FALSE
	
	return(object)
}


transformKernel <- function(kern, linkinfo, type=c("none", "simple", "extended"), linkfun=function(x) x^(1/6), 
							notlinkfun=function(x) {1 - (1-x)^(1/6)}) {
	# transform kernel using linkinfo object.
	# "simple" : transform only similarities concerned by linking functions
	# "extended" : attach all elements that are not in the U and L relations to their closest neighbour in the relations (1NN)
	# 			and augment the relation with the set of elements associated to this closest neighbour, and link to that 
	#			closest element -> just add a "link" with the closest element : transtivity does the rest of the job.
	# linking functions may be parametrized away from defaults.
	type <- match.arg(type)
	relatedElts <- numeric(0)
	unrelatedElts <- numeric(0)
	
	# for debug purposes
	#savelink <- linkinfo

	indexes <- function(vec) {
		vec <- which(vec==TRUE)
		relatedElts <<- c(relatedElts, vec) # modification in the enclosing lexical scope, ie calling environment
		return(vec)
	}
	
	maxRelated <- function(ind) {
		maxind <- which.max(kern[ind,relatedElts])
		maxind <- relatedElts[maxind] # cast back to original index
		linkinfo <<- add(linkinfo, ind, maxind, type="link") # direct modification of linkinfo
		return(maxind)
	}

	
	# specific if "extended" : find those elements which are not in any relation,
	# and "attach" them to their closest element that is in a relation.
	if (type == "extended") {
		diag(kern) <- 0 # computation based on max similarity : exclude self.
		apply(linkinfo$link, 1, indexes)
		apply(linkinfo$notlink, 1, indexes)
		relatedElts <- unique(relatedElts)
		unrelatedElts <- 1:length(linkinfo)
		if (length(relatedElts > 0)) {
			unrelatedElts <- unrelatedElts[-relatedElts]
			apply(as.array(unrelatedElts), 1, maxRelated)
		}
		diag(kern) <- 1 # reestablish correct kernel
	}
	
	if (!check(linkinfo)) {
		#browser() # debug purposes
		stop("linkinfo is not consistent")
	}
	
	# apply transitivity relations
	templink <- transform(linkinfo)
	
	#browser()
	
	# transform only if appropriate
	if (type != "none") {
		kern[templink$link] <- linkfun(kern[templink$link])
		kern[templink$notlink] <- notlinkfun(kern[templink$notlink])
	}
	return(kern)
}




distortions <- function(origdata, projdata) {
	# standardize both data sets
	Reg <- 10^(-3) * dim(origdata)[1]
	origmeans <- apply(origdata, 2, mean)
	projmeans <- apply(projdata, 2, mean)
	origsds <- apply(origdata, 2, sd)
	projsds <- apply(projdata, 2, sd)
	origsds[origsds==0] <- 1 # manage NULL case
	projsds[projsds==0] <- 1
	origdata <- sweep(origdata, 2, origmeans, "-")
	projdata <- sweep(projdata, 2, projmeans, "-")
	origdata <- sweep(origdata, 2, origsds, "/")
	projdata <- sweep(projdata, 2, projsds, "/")
	
	# compute distance matrices
	origmat <- as.matrix(dist(origdata))
	projmat <- as.matrix(dist(projdata))
	# normalize in [0,1]
	origmat <- origmat / max(origmat)
	projmat <- projmat / max(projmat)
	
	# compute D^+ and D^- matrices
	posdiffs <- negdiffs <- diffs <- origmat - projmat
	posdiffs[posdiffs<0] <- 0
	negdiffs[negdiffs>0] <- 0
	negdiffs <- -negdiffs
	
	# mu aggregates, and results
	mucompress <- apply(posdiffs, 1, sum) # aggregate all columns for each row-element
	mustretch <- apply(negdiffs, 1, sum)
	
	# if all distortions are almost 0, then we may have significantly \neq 0 normalized dists -> problem.
	# see with michael's paper. maybe normalize with some maximal/minimal pairwise distance.
	# propose to "regularize" it - we use a Laplace regularization. In the context of summed pairwise distances, n.10^(-3) seems a reasonable value.
	rescompress <- apply(as.data.frame(mucompress), 1, 
							function(x) (x - min(mucompress)) / (max(mucompress) - min(mucompress) + Reg))
	resstretch <- apply(as.data.frame(mustretch), 1, 
							function(x) (x - min(mustretch)) / (max(mustretch) - min(mustretch)+ Reg))
	return(list(compress=rescompress, stretch=resstretch))
}








