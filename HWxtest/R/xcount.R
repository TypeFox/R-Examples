# xcount
# (c) William R. Engels, 2014
# 'Exact Tests For Hardy-Weinberg Proportions', 2009, Genetics, 183, pp1431-1441


#' Find Exact Number of Genotype Tables
#' 
#' Use \code{xcount} to determine the exact number of tables (i.e., genotype numbers) for a given set of allele counts. This method enumerates all tables, and is best when the total number is less than 10^10 or so. This function is mostly called by \code{\link{hwx.test}} rather than directly by the user. If the number of tables is too large to enumerate with this method, use \code{\link{acount}} for an approximation.
#' 
#' @param m vector containing the numbers of alleles of each type. Length must be at least 2 and all must be non-negative integers. It can also be a \code{matrix} of genotype counts.
#' @param safety Stop execution if the approximate table number obtained from \code{\link{acount}} is more than this cutoff.
#' @param safeSecs Time limit in seconds. Another safety feature to prevent getting stuck in a too-long computation
#' 
#' @return The exact number of tables
#' 
#' @examples
#' # Allele counts from human Rh locus. Guo and Thompson, 1992, Figure 1
#' #
#' alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
#' xcount(alleles)
#' 
#' @references The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
#' 
#' @seealso \code{\link{hwx.test}}, \code{\link{acount}}


# dyn.load("~/DropBox/HWxtest/pkg/src/HWxcount.so")


#' @useDynLib HWxtest
#' @export
xcount <- 
function(m, safety = 1e10, safeSecs = 10) {
	UseMethod("xcount")
}


#' @export
xcount.integer <- 
function(m, safety = 1e10, safeSecs = 10) {
	m <- m[m>0]
	if(length(m) < 2) return(1);
	if(any(m < 0)) stop("\nAllele counts must be nonnegative\n");
	ac <- acount(m);
	if(ac > safety) stop("\nToo many to count. Try increasing safety paramater or call acount\n");
	xc <- -1;
	value <- .C("xcount",
		counts=as.integer(sort(m, decreasing=T)),
		nAlleles = as.integer(length(m)),
		tableCount=as.double(xc),
		safeSecs=as.integer(safeSecs)
		,PACKAGE="HWxtest"
		);
	n <- value$tableCount;
	if(n < 0) {
		warning("\nOperation timed out after ", -n, 
		"tables.\nThe true count is greater than that.\nNegative value for count reported.\nTime limit is defined by parameter safeSecs.")
	}
	value$tableCount;
}

#' @export
xcount.matrix <- 
function(m, safety = 1e10, safeSecs = 10) {
	m <- alleleCounts(m);
	xcount(m, safety=safety, safeSecs=safeSecs)
}


#' @export
xcount.table <- 
function(m, safety = 1e10, safeSecs = 10) {
	m <- unclass(m);
	xcount(m, safety=safety, safeSecs=safeSecs)
}

#' @export
xcount.numeric <- 
function(m, safety = 1e10, safeSecs = 10) {
	m <- as.integer(m);
	xcount(m, safety=safety, safeSecs=safeSecs)
}

#' @export
xcount.genotype <- 
function(m, safety = 1e10, safeSecs = 10) {
	requireNamespace("genetics")
	tab <- table(factor(genetics::allele(m, 1), levels = genetics::allele.names(m)), factor(genetics::allele(m, 2), levels = genetics::allele.names(m)));
	m <- alleleCounts(unclass(t(tab)));
	xcount(m, safety=safety, safeSecs=safeSecs)
}


#' @export
xcount.logical <- 
function(m, safety = 1e10, safeSecs = 10) {return(NA)}
