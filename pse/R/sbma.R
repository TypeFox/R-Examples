#' Calculates the Symetrized Blest Measure of Agreement between two samples
#'
#' The Symetrized Blest Measure of Agreement is an alternative measure of
#' rank correlation (similar to Kendall's Tau and Spearman's Rho). This correlation
#' measure is more sensitive to changes in the order of the first elements of a vector 
#' (see examples).
#'
#' This	function calculates the SBMA between two samples or two \code{\link{LHS}} objects.
#' In the second case, what is compared is the values of the "prcc" component of
#' each Hypercube.
#'
#' @param sample1 The first vector or LHS object to be compared.
#' @param sample2 The second vector or LHS object to be compared.
#' @param absolute Logical. Should the absolute values of sample1 and sample2 be used 
#' in the calculation?
#' @param \dots Additional arguments. 
#' @references
#'  Maturi, T.A. and Elsayigh, A. 2010. A comparison of correlation coefficients 
#'  via a three-step bootstrap approach. \emph{Journal of Mathematics Research} 2(2): 3-10.
#' @examples
#' # SBMA is only affected by the rank of the values inside each vector
#' sbma(c(1,2,3,4), c(2,3,4,5))
#' # Changes in the first positions: high impact on the SBMA
#' sbma(c(1,2,3,4), c(2,1,3,4))
#' cor(c(1,2,3,4), c(2,1,3,4), method="spearman")
#' # Changes in the last positions: low impact on the SBMA
#' sbma(c(1,2,3,4), c(1,2,4,3))
#' cor(c(1,2,3,4), c(1,2,4,3), method="spearman")
#' @export
sbma <- function (sample1, sample2, absolute=TRUE, ...) UseMethod("sbma")

#' @export
#' @rdname sbma
sbma.default <- function (sample1, sample2, absolute=TRUE, ...) {
	x1 <- sample1; x2<-sample2;
	if (absolute) {x1 <- abs(x1); x2 <- abs(x2);}
	if (length(x1) != length(x2))
		stop("Sample sizes must be the same!");
	n <- length(x1);
	R <- rank(x1);
	S <- rank(x2);
	v1 <- -(4*n+5)/(n-1);
	v2 <- 6/(n^3-n);
	sbc <- v1+v2*(sum(R*S*(4-(R+S)/(n+1))));
	return (sbc)
}

#' @export
#' @rdname sbma
sbma.LHS <- function(sample1, sample2, absolute, ...) {
	sb <- array();
	for (i in 1:dim(get.results(sample1))[2])
		sb[i] <- sbma(sample1$prcc[[i]]$PRCC$original, sample2$prcc[[i]]$PRCC$original);
	return (sb);
}   
