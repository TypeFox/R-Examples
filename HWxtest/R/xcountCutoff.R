# xcountCutoff
# (c) William R. Engels, 2014
# 'Exact Tests For Hardy-Weinberg Proportions', 2009, Genetics, 183, pp1431-1441


#' Determine immediately whether number of tables is over a limit
#' 
#' Calling \code{scountCutoff} gives you a quick answer to whether the number of tables is over a given cutoff. It is useful in deciding whether to analyze a data set with \code{\link{xtest}} or \code{\link{mtest}}. This function is used by \code{\link{hwx.test}} and not normally called directly by the user.
#' 
#' 
#' @param m vector containing the numbers of alleles of each type. It can also be a matrix of genotype counts, but not a vector of genotype counts.
#' @param cutoff Is the number of tables above or below this value?
#' 
#' 
#' @return TRUE or FALSE depending on whether the table count is above or below \code{cutoff}
#' 
#' @examples
#' #
#' alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
#' if(xcountCutoff(alleles)) cat("There are too many tables")


#' @useDynLib HWxtest

#' @export
xcountCutoff <- 
function(m, cutoff=1e7) {
	UseMethod("xcountCutoff")
}

#' @export
xcountCutoff.integer <- 
function(m, cutoff=1e7) {
	m <- m[m!=0]
	if(length(m) < 2) return(TRUE);
	if(any(m < 0)) stop("\nAllele counts must be nonnegative\n");
	if(acount(m) > 1e11) return(TRUE);
	value <- .C("xcount",
		counts=as.integer(sort(m, decreasing=T)),
		nAlleles = as.integer(length(m)),
		tableCount=as.double(cutoff),
		safeSecs=as.integer(5)
		,PACKAGE="HWxtest"
	);
		n <- value$tableCount;
		if(n < 0) return(TRUE);
		return(FALSE)
}

#' @export
xcountCutoff.matrix <- 
function(m, cutoff=1e7) {
	m <- alleleCounts(m);
	xcountCutoff(m, cutoff=cutoff)
}

#' @export
xcountCutoff.numeric <-
function(m, cutoff=1e7) {
	m <- as.integer(m);
	xcountCutoff(m, cutoff=cutoff)
}