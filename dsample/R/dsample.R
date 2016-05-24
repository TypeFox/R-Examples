#' @rdname dsample
#' @title Random Samples Generation Through The Wang-Lee and Fu-Wang Algorithms
#' @description \code{sample.wl} generates a sample of specified size \code{n} from the target density funciton (up to a normalizing constant) based on the Wang-Lee algorithm
#' @param expr expression
#' @param rpmat matrix containing random points for discretization
#' @param  nk     a positive integer, the number of contours.  See \sQuote{Details}.
#' @param  n      a non-negative integer, the desired sample size. 
#' @param  wconst a real number between 0 and 1.  See \sQuote{Details}.
#' @details  \code{X} has the number of rows equals to the number of discrete base points. In each row, the first element contians the funcitonal value of the target density and the rest elements are the coordinates at which the density is evaluated.  
#' \code{wconst} is a constant for adjusting the volumn of the last contour.
#' @return \code{sample.wl} gives the drawn sample as a \code{data.frame} with number of rows equals the specified size \code{n} and number of columns equals \code{ncol(x)-1}.
#' @references
#' Wang, L. and Lee, C.H. (2014). Discretization-based direct random sample generation. Computational Statistics and Data Analysis, 71, 1001-1010. 
#' Lee, C.H. (2009). Efficient Monte Carlo Random Sample Generation through Discretization, MSc thesis, Department of Satistics, University of Manitoba, Canada
#' Wang, L. and Fu, J. (2007). A practical sampling approach for a bayesian mixture model with unknown number of components. Statistical Papers, 48(4):631-653.
#' Fu, J. C. and Wang, L. (2002). A random-discretization based Monte Carlo sampling method and its application. Methodology and Computing in Applied Probability, 4, 5-25.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}, Liqun Wang \email{liqun.wang@@umanitoba.ca}
#' @keywords sampling, discretization
#' @examples 
#' ## The following example is taken from West (1993, page 414).
#' ## West, M. (1993). Approximating posterior distributions by mixture.
#' ##   Journal of the Royal Statistical Society - B, 55, 409-422.
#' 
#' expr <- expression((x1*(1-x2))^5 * (x2*(1-x1))^3 * (1-x1*(1-x2)-x2*(1-x1))^37)
#' sets <- list(x1=runif(1e5), x2=runif(1e5))
#' smp <- dsample(expr=expr, rpmat=sets, nk=1e4, n=1e3)
#' 
#' ##
#' ## More accurate results can be achieved by increasing the number 
#' ## of dicretization points and the number of contours.  
#' @export 
dsample <- function(expr, rpmat, n=1e3, nk=1e4, wconst){
	
	# option 1. random discretization (currently used)
	# option 2. grid discretization (approximation)
	
	variables <- all.vars(expr)
	if(missing(wconst)) wconst <- 1
	cl <- match.call(expand.dots=FALSE)
	stopifnot(is.expression(expr), is.list(rpmat), is.numeric(nk), is.numeric(n), is.numeric(wconst))	

	y <- eval(expr=expr, envir=rpmat)
	fmla <- stats::as.formula(paste("y~", expr, sep=""))
	
	X <- as.data.frame(rpmat)
	yX <- cbind(y, X)
	yX <- yX[which(y>0),] # data.frame
	yX <- yX[order(yX$y, decreasing=TRUE), ]
	
	cnt <- graphics::hist(yX$y, breaks=seq(from=min(yX$y), to=max(yX$y), length.out=nk+1), plot=FALSE)
	cnames <- paste("e", seq_len(nk), sep="") # contour name (in order)
	yX$cid <- cut(yX$y, cnt$breaks, include.lowest=TRUE)
	yX$cnt.name <- cnames[match(yX$cid, names(table(yX$cid)))]
	cnt.counts <- cnt$counts
	cnt.mids <- cnt$mids
	names(cnt.mids) <- names(cnt.counts) <- cnames
	
	gpdf <- rev(cnt.counts * cnt.mids)
	gpdf[nk] <- wconst*gpdf[nk]
	rev.cntc <- rev(cnt.counts)
	
	cdf <- cumsum(gpdf)/sum(gpdf)
	cumcdf <- cumsum(cdf)
	
	# counting how many samples are needed from each contour
	pptns <- graphics::hist(stats::runif(n), breaks=c(0,cdf), plot=FALSE)$counts
	names(pptns) <- rev(paste("e", seq_len(nk), sep=""))
	
	# sampling from each contour 
	scnt <- mapply(FUN=sample, MoreArgs=list(replace=TRUE), rev.cntc, pptns)
	idx <- unlist( mapply("+", as.list( c(0, cumsum(rev.cntc))[-(nk+1)] ), scnt) )
	yX <- yX[order(yX$y, decreasing=TRUE)[idx], ]
	
	robj <- list(formula=fmla, expr=expr, yX=yX, X=yX[all.vars(expr)], cnt.counts=cnt.counts, cnt.mids=cnt.mids, gpdf=gpdf, cdf=cdf, cumcdf=cumcdf, pptns=pptns, scnt=scnt, idx=idx)
	class(robj) <- "dsample"
	return(robj)
}


#' @rdname summary_dsample
#' @title Generating Basic Summary Statistics of Marginal Distributions
#' @description  Producing basic summary statistics (the mean, the standard deviation and the first five modes) from the sample drawn via either the Fu-Wang algorithm or the Wang-Lee algorithm, for all marginal distributions of the target distribution.
#' @param object a \code{data.frame}, contains the sample drawn via either the Fu-Wang algorithm or the Wang-Lee algorithm 
#' @param n the first n samples
#' @param ... more arguments
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}, Liqun Wang \email{liqun.wang@@umanitoba.ca}
#' @export 
summary.dsample <- function(object, n=5, ...) {

	stopifnot(inherits(object, "dsample"))
	fmla <- object$formula
	
	X <- object$X
	hc <- stats::hclust(stats::dist(X))
	grp <- stats::cutree(hc, ...)
	cdf <- object$cdf
	
	means <- colMeans(X)
	stdevs <- do.call(c, lapply(X, stats::sd, na.rm=TRUE))
	modes <- cbind(X)[1:n,]
	
	robj <- list(formula=fmla, means=means, stdevs=stdevs, modes=modes, hc=hc, grp=grp, X=X, cdf=cdf)
	class(robj) <- "dsample"
	return(robj)
}


#' @rdname plot
#' @title Plot dsample objects
#' @description Plot dsample objects 
#' @param x dsample object.
#' @param ... arguments passing functions inside.
#' @export 
plot.dsample <- function(x, ...){

	X <- x$X
	cdf <- x$cdf
	grp <- x$grp

	graphics::par(mfrow=c(2,2))
	graphics::plot(cdf, main="CDF", xlab="E", ylab="F(E)", cex=0.5)
	graphics::plot(X, cex=0.5, main="Scatter and Contour Plots", xlab="x1", ylab="x2", col=grp)
	density <- MASS::kde2d(X[,1], X[,2], n=1e3)
	graphics::contour(density, nlevels=5, add=TRUE)
	graphics::hist(X[,1], main="Histogram", ylab="density", xlab=expression(x[1]), prob=TRUE, breaks=20)
	graphics::hist(X[,2], main="Histogram", ylab="density", xlab=expression(x[2]), prob=TRUE, breaks=20)
}
