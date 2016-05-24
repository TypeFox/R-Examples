.sort.missing <- function(x, x_nonmiss){
	## Sort the observations based on their missingness for faster computation (in partial mahalanobis distance)
	## Added Mar 17, 2012
	miss.group <- apply( x_nonmiss, 1, paste, collapse='')
	id.order <- order(miss.group, decreasing=TRUE)
	miss.group.unique <- do.call(rbind, strsplit(sort(unique(miss.group), decreasing=TRUE),""))
	miss.group.unique <- t(apply(miss.group.unique, 1, as.numeric))
	miss.group.counts <- as.numeric(table(miss.group)[unique(miss.group[ id.order])])
	miss.group.n <- length(miss.group.counts)
	miss.group.p <- rowSums(miss.group.unique)
	p <- ncol(x); n <- nrow(x)
	miss.group.obs.col <- sapply( 1:miss.group.n, function(i) {
		tmp <- rep(0,p)
		if( miss.group.p[i] > 0 ) tmp[1:miss.group.p[i]] <- which(miss.group.unique[i,]==1)
		tmp
	})
	miss.group.mis.col <- sapply( 1:miss.group.n, function(i) {
		tmp <- rep(0,p)
		if( p-miss.group.p[i]> 0)tmp[1:(p-miss.group.p[i])] <- which(miss.group.unique[i,]==0)
		tmp
	})

	miss.group.obs.col <- t(miss.group.obs.col)
	miss.group.mis.col <- t(miss.group.mis.col)
	
	x.miss.group.match <- rep(1:miss.group.n, miss.group.counts)

	list(id.order=id.order, 
		id.ro=order(id.order),
		x=x[id.order,],
		x_nonmiss=x_nonmiss[id.order,],
		pu=rowSums(x_nonmiss[id.order,]),
		x.miss.group.match=x.miss.group.match,
		miss.group.unique=miss.group.unique, 
		miss.group.counts=miss.group.counts, 
		miss.group.obs.col=miss.group.obs.col,
		miss.group.mis.col=miss.group.mis.col,
		miss.group.p=miss.group.p,
		miss.group.n=miss.group.n) 
}


## Compute partial mahalanobis distance 
## this is for the case when the data set is not sorted 
partial.mahalanobis <- function(x, mu, Sigma){
	xcall <- match.call()

	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")

	## drop all rows with missing values (!!) :
	n <- nrow(x); p <- ncol(x) 
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")
		
	## Rows with at least one observed
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	x.orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]	
	
	## get missing pattern and sorted matrix
	x_sort <- .sort.missing(x, x_nonmiss)
	
	## Compute pmd
	if( any( x_nonmiss == 0 ) ){
		x_cent <- sweep(x_sort$x, 2, mu, "-")
		pmd.tmp <- .partial.mahalanobis.Rcpp( x_cent, Sigma, x_sort$miss.group.unique, x_sort$miss.group.counts)
		pmd.tmp <- pmd.tmp[ x_sort$id.ro ]
	} else{
		pmd.tmp <- mahalanobis( x_sort$x, mu, Sigma)
		pmd.tmp <- pmd.tmp[ x_sort$id.ro ]
	}

	pmd <- rep(NA, nrow(x.orig))
	pmd[ok] <- pmd.tmp
	
	pu <- rowSums( !is.na(x.orig))
	pmd.adj <- qchisq( pchisq( pmd, df=pu, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pu == p)] <- pmd[ which(pu==p) ]	
	
	new("CovRobMiss", 
		mu = mu,
		S = Sigma,
		call = xcall,
		estimator = "unknown",
		x = x.orig, 
		pmd = pmd,
		pmd.adj = pmd.adj,
		p=p,
		pu = pu)	
}

