###################################################################
## Rcpp version of log-density of multivariate normal
.fast.ldmvnorm.Rcpp <- function(x_mu_diff, Sigma, miss_group_unique, miss_group_counts){
	res <- tryCatch( .Call("fast_mvnorm_density", x_mu_diff, Sigma, miss_group_unique, miss_group_counts),
		"std::range_error" = function(e){
		conditionMessage( e ) } )
	if( is.character( res ) ) stop(res)
	return( c(res) )
}

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



ldmvnorm <- function(x, mu, Sigma, onNA=0){
	if( missing(x) ) stop("'x' missing")
	xcall <- match.call()

	## Check x
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")

	## Check mu and Sigma
	n <- nrow(x); p <- ncol(x) 
	if (missing(mu)) mu <- rep(0, p)
    if (missing(Sigma)) Sigma <- diag(p)
	if( length(mu) != p) stop("'x' and 'mu' have non-conforming size")
	if( ncol(Sigma) != p) stop("'x' and 'Sigma' have non-conforming size")

	## drop all rows with missing values (!!) :
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")
	if( all(pp == p) ) return( dmvnorm(x=x, mean=mu, sigma=Sigma, log=TRUE) )
	
	## Rows with at least one observed
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	x.orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]	
	
	## get missing pattern and sorted matrix
	x_sort <- .sort.missing(x, x_nonmiss)
	
	## Compute log density
	x_cent <- sweep(x_sort$x, 2, mu, "-")
	dens.tmp <- .fast.ldmvnorm.Rcpp( x_cent, Sigma, x_sort$miss.group.unique, x_sort$miss.group.counts)
	dens.tmp <- dens.tmp[ x_sort$id.ro ]
	dens <- rep(NA, nrow(x.orig))
	dens[ok] <- dens.tmp

        ## onNA 

        dens[is.na(dens)]=onNA 
	return(dens)
}

