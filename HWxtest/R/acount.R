# acount
# (c) William R. Engels, 2014
# 'Exact Tests For Hardy-Weinberg Proportions', 2009, Genetics, 183, pp1431-1441


#' Find Approximate Number of Genotype Tables
#' 
#' Use \code{acount} to obtain the approximate number of genotype tables for a given set of allele counts. This method uses a normal approximation and is much faster than enumerating the tables with \code{\link{xcount}} but not as accurate.
#' 
#' @param m vector containing the numbers of alleles of each type. Length must be at least 2. All items are positive integers. It can also be a matrix of genotype counts, an object of type \code{genotype}, but not a vector of genotype counts.
#' 
#' @return The approximate number of tables.
#' 
#' @examples
#' # Allele counts from human Rh locus. Guo and Thompson, 1992, Figure 1
#' #
#' alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
#' acount(alleles)
#' # This approximation may be compared with the exact value of 250552020
#' #
#' ld <- c(6329, 319, 47, 2773, 75, 6702, 14, 2, 333)
#' acount(ld)
#' #
#' # This is an example where the number of tables is too large for a full enumeration.
#' 
#' @references The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
#' 
#' @seealso \code{\link{hwx.test}}, \code{\link{xcount}}



#' @export acount
acount <- 
function(m){
	UseMethod("acount")
}


#' @export
acount.integer <- 
function(m) {
	m <- m[m != 0]
	if(length(m) < 2) return(1);
	if(any(m < 0)) stop("\nAllele counts must be nonnegative\n");
	n <- sum(m)/2;
	k <- length(m);
	summ2 <- sum(m^2);
	b <- k * (k+1)/2 -1;
	va <- n * b * (n+b+1)/((b+1)*(b+1)*(b+2));
	vm <- (k+1) * va;
	q <- ((k-1)/(vm * k)) * (summ2 - (4 * n * n/k));
	lnPm <- log(sqrt(k))  +  ((k-1)/2) * log((k-1)/(2*pi*k*vm)) - q/2;
	lns <- lchoose(n + b, b);
	exp(lns + lnPm)
}


#' @export
acount.matrix <- 
function(m) {
	m <- alleleCounts(m);
	acount.integer(m)
}

#' @export
acount.table <- 
function(m) {
	acount(unclass(m))
}

#' @export
acount.numeric <- 
function(m) {
	m <- as.integer(m);
	acount.integer(m)
}


#' @export
acount.genotype <- 
function(m) {
	requireNamespace("genetics")
	tab <- table(factor(genetics::allele(m, 1), levels = genetics::allele.names(m)), factor(genetics::allele(m, 2), levels = genetics::allele.names(m)));
	m <- alleleCounts(unclass(t(tab)));
	acount.integer(m)
}


#' @export
acount.logical <- 
function(m) {return(NA)}