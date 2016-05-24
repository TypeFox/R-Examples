# Full enumeration test for HW
# (c) William R. Engels, 2014


#' Performs an exact test with full enumeration for Hardy-Weinberg proportions.
#' 
#' Given a set of genotype counts, \code{xtest} examines all possible outcomes with the same set of allele counts. For each table, it computes four test statistics and compares them with the observed values. It returns the total probability of all tables with test statistics as \dQuote{extreme} or more so than the observed. It can also plot a histogram of one of the statitistics if \code{histobins} is greater than zero. More about these four test statistics and other information can be found in the vignette. This function will not normally be called directly. Instead, \code{\link{hwx.test}} calls either \code{xtest} or \code{\link{mtest}} depending on which method is to be used.
#' 
#' @param c A matrix containing the genotype counts. It should be a square matrix, but only the lower-left half is used.
#' @param statName can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.
#' @param histobins If 0 no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If set to a number greater than 1, a histogram with \code{histobins} is plotted.
#' @param histobounds A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{xtest} will compute reasonable bounds to include most of the distribution.
#' @param showCurve whether to show a blue curve indicating the asymptotic (chi squared) distribution. This only works for \code{LLR} and \code{Chisq}
#' @param safeSecs After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.
#' @param detail Determines how much detail is printed. If set to 0, nothing is printed (useful if you use \code{xtest} programmatically.).
#' 
#' @return \code{xtest} returns a list components
#' \item{$ Pvalues}{The four computed P values corresponding to the test statistics: \code{LLR}, \code{Prob}, \code{U} and \code{Chisq} in that order.}
#' \item{$ observed}{The four observed statistics in the same order as above}
#' \item{$ tableCount}{The number of tables examined during the calculation}
#' \item{$ ntrials}{placeholder}
#' \item{$ genotypes}{The input matrix of genotype counts}
#' \item{$ alleles}{The allele counts \eqn{m} corresponding to the input genotype counts}
#' \item{$ statID}{Which test statistic was used if a histogram was plotted}
#' \item{$ histobins}{If greater than zero, the number of bins to use for the histogram}
#' \item{$ histobounds}{The lower and upper limits of the test statistic in the histogram}
#' \item{$ histoData}{Vector of \eqn{histobins} values for the histogram}
#' \item{$ showCurve}{Whether the asymptotic curve should be plotted with the histogram}

#' 
#' @references The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
#' 
#' @seealso \code{\link{hwx.test}}


#' @useDynLib HWxtest

xtest <- 
function(c, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	statNames <- c("LLR", "Prob", "U", "Chisq");
	statID <- which(statNames==statName);
	m <- alleleCounts(c); names(m) <- rownames(c);
	if(histobins == 1){ histobins  <- 500};   #The default is 500 bins
	ostats <- c(observedLLR(c), observedProb(c), observedU(c), observedX2(c));
	names(ostats) <- statNames;
	if(histobounds[1]==histobounds[2] && histobins) histobounds <- defaultHistobounds(ostats, statID, m);
	x <- .C("xtest",
			m=as.integer(sort(m, decreasing=T)),
			nAlleles=as.integer(length(m)),
			observed=as.double(ostats),
			Pvalues=as.double(double(4)),
			statID=as.integer(statID-1),
			histobins=as.integer(histobins),
			histobounds=as.double(histobounds),
			histoData=as.double(double(histobins + 1)),
			safeSecs=as.integer(safeSecs),
			tableCount=as.double(0)
			,PACKAGE="HWxtest"
			);
	if(histobins) { # Convert histo counts to densities
		binwidth <- (histobounds[[2]] - histobounds[[1]])/histobins
		x$histoData  <- x$histoData/binwidth
	}		

	names(x$Pvalues) <- statNames;
	se <- rep(NA, 4); names(se) <- statNames;
	return = list(Pvalues=x$Pvalues, 
		observed=ostats, 
		tableCount=x$tableCount, 
		ntrials= NA, 
		genotypes=c, 
		alleles=m, 
		SE=se,
		statID=statID,
		histobins=histobins,
		histobounds=histobounds,
		histoData=x$histoData,
		showCurve=showCurve)

}