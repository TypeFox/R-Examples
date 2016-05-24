# Compute statistics for observed table in HW test
# (c) William R. Engels, 2014


#' Compute observed statistics for a genotype count matrix
#' 
#' #' Four measures of fit to Hardy-Weinberg for a given set of genotype counts may be computed.
#' \itemize{
	#' \item \code{observedProb} The probability of the observed set under the HW null and with the allele counts fixed.
	#' \item \code{observedLLR} The log-likelihood ratio of the observed set
	#' \item \code{observedU} The observed U-score. Positive values indicate an excess of homozygotes and negative ones imply too many heterozygotes
	#' \item \code{observedX2} The classical \dQuote{chi-squared} statistic
#' }
#' 
#' 
#' @param c Matrix of observed genotype counts. Each number should be a non-negative integer, and matrix is \eqn{k x k}.
#' @param returnExpected Used in \code{observedX2} to indicate whether a matrix of expected numbers should be returned instead.
#' 
#' @return the observed statistic
#' 
#' @examples
#' t <- vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
#' observedStats <- c(observedProb(t), observedLLR(t), observedU(t), observedX2(t))


#' @rdname observedStatistics
#' @export
observedProb <- 
function(c){
	c <- fillUpper(c);
	m <- alleleCounts(c);
	d <- sum(diag(c));
	n <- sum(m)/2;
	k <- n * log(2) + lgamma(n+1)- lgamma(2 * n + 1) + sum(lgamma(m+1));
	a <- matrix.to.vec(c);
	p <- exp(k - sum(lgamma(a+1)) - d*log(2));
}
#' @rdname observedStatistics
#' @export
observedLLR <- 
function(c){
	c <- fillUpper(c);
	m <- alleleCounts(c);
	d <- sum(diag(c));
	n <- sum(m)/2;
	a <- matrix.to.vec(c);
	k <- sum(m * log(m), na.rm=T) - (n+d)*log(2) - n*log(n);
	llr <- k - sum(a * log(a), na.rm=T);
}

#' @rdname observedStatistics
#' @export
observedU <- 
function(c){
	c <- fillUpper(c);
	m <- alleleCounts(c);
	hz <- diag(c);
	n <- sum(m)/2;
	u <- n*(2*sum(hz/m) - 1);	
}

#' @rdname observedStatistics
#' @export
observedX2 <- 
function(c, returnExpected=F){
	c <- fillUpper(c);
	k <- dim(c)[1];
	m <- alleleCounts(c);
	d <- sum(diag(c));
	n <- sum(m)/2;
	a <- matrix.to.vec(c);
	ae <- matrix(NA,k,k)
	for(i in 1:k) for(j in 1:k) ae[i,j] <- m[i]*m[j]/(2*n);
	for(i in 1:k)ae[i,i] <- ae[i,i]/2;
	ve <- matrix.to.vec(ae);
	if(returnExpected) return(vec.to.matrix(ve));
	sum((a-ve)^2/ve)
}

# Find the probability of a perfectly-fitting set of genotype counts for the given allele counts.
perfectProb <- 
function(c){
	m <- alleleCounts(c);
	k <- dim(c)[1];
	n <- sum(m)/2;
	ae <- matrix(0,k,k);
	for(i in 1:k) for(j in 1:k) ae[i,j] <- m[i]*m[j]/(2*n);
	for(i in 1:k)ae[i,i] <- ae[i,i]/2;
	#ae <- round(ae);
	observedProb(ae);
	ae
}