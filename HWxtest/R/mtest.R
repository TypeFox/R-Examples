# Monte Carlo test for HW
# (c) William R. Engels, 2014


#' Performs an \dQuote{exact} test using Monte Carlo trials for Hardy-Weinberg proportions
#' 
#' Given a set of genotype counts, \code{mtest} examines a large number of possible outcomes with the same set of allele counts. For each table, it computes four test statistics and compares them with the observed values. It returns the total probability of all tables with test statistics as \dQuote{extreme} or more so than the observed. It can also plot a histogram of one of the statitistics if \code{histobins} is greater than zero. More about these four test statistics and other information can be found in the vignette. This function will not usually be called directly by the user. Instead, call \code{\link{hwx.test}} with \code{method} set to either \dQuote{auto} or \dQuote{monte}.
#' 
#' @param c A matrix containing the genotype counts. It should be a square matrix, but only the lower-left half is used.
#' @param ntrials the number of random trials to perform
#' @param statName can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.
#' @param histobins If 0 no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If set to a number greater than 1, a histogram with \code{histobins} is plotted.
#' @param histobounds A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{mtest} will compute reasonable bounds to include most of the distribution.
#' @param showCurve whether to show a blue curve indicating the asymptotic (chi squared) distribution. This only works for \code{LLR} and \code{Chisq}
#' @param safeSecs After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.
#' @param detail Determines how much detail is printed. If it is set to 0, nothing is printed (useful if you use \code{mtest} programmatically.).
#' 
#' @return \code{mtest} returns a list components
#' \item{$ Pvalues}{The four computed P values corresponding to the test statistics: \code{LLR}, \code{Prob}, \code{U} and \code{Chisq} in that order.}
#' \item{$ tableCount}{placeholder}
#' \item{$ SE}{Standard errors for the P values. These come from the binomial.}
#' \item{$ observed}{The four observed statistics in the same order as above}
#' \item{$ ntrials}{The number of tables examined during the calculation}
#' \item{$ genotypes}{The input matrix of genotype counts}
#' \item{$ alleles}{The allele counts \eqn{m} corresponding to the input genotype counts}
#' \item{$ statID}{Which test statistic was used if a histogram was plotted}
#' \item{$ histobins}{If greater than zero, the number of bins to use for the histogram}
#' \item{$ histobounds}{The lower and upper limits of the test statistic in the histogram}
#' \item{$ histoData}{Vector of \eqn{histobins} values for the histogram}
#' \item{$ showCurve}{Whether the asymptotic curve should be plotted with the histogram}
#' 
#' 
#' @references The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
#' 
#' @seealso \code{\link{hwx.test}}
#' 


#' @useDynLib HWxtest

mtest <- 
function(c, ntrials=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	statNames <- c("LLR", "Prob", "U", "Chisq");
	statID <- which(statNames==statName);
	m <- alleleCounts(c); names(m) <- colnames(c)
	if(histobins == 1){ histobins  <- 500};   #The default is 500 bins
	ostats <- c(observedLLR(c), observedProb(c), observedU(c), observedX2(c));
	names(ostats) <- statNames;
	if(histobounds[1]==histobounds[2] && histobins) histobounds <- defaultHistobounds(ostats, statID, m);
	runif(20); # Just to get random numbers started
	x <- .C("mtest",
			m=as.integer(sort(m, decreasing=T)),
			nAlleles=as.integer(length(m)),
			observed=as.double(ostats),
			Pvalues=as.double(double(4)),
			statID=as.integer(statID-1),
			histobins=as.integer(histobins),
			histobounds=as.double(histobounds),
			histoData=as.double(double(histobins + 1)),
			safeSecs=as.integer(safeSecs),
			ntrials=as.double(ntrials)
			,PACKAGE="HWxtest"
			);
	if(histobins) { # Convert histo counts to densities
		binwidth <- (histobounds[[2]] - histobounds[[1]])/histobins
		x$histoData  <- x$histoData/(ntrials * binwidth)
	}		
	names(x$Pvalues) <- statNames;
	se <- sqrt(abs(x$Pvalues * (1-x$Pvalues)/x$ntrials))
	names(se) <- statNames
	return = list(Pvalues=x$Pvalues, 
		observed=ostats, 
		tableCount=NA, 
		ntrials= x$ntrials, 
		genotypes=c, 
		alleles=m, 
		SE=se,
		statID=statID,
		histobins=histobins,
		histobounds=histobounds,
		histoData=x$histoData,
		showCurve=showCurve)

}