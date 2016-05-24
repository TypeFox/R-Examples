# hwx.test  Generic function to test for HW by either full enumeration or by Monte Carlo
# (c) William R. Engels, 2014
# Accepts matrix, vector, table, genotype, etc. 


#' Test for HW by either full enumeration or Monte Carlo.
#' 
#' The \code{hwx.test()} function is the main function of the \code{HWxtest} package. This function produces a valid test for Hardy-Weinberg frequencies for virtually any set of genotype counts. It will use either a full-enumeration method in which all possible tables with the same allele numbers are examined, or a Monte Carlo test where a large number of random tables is examined. To decide which to use, it calls \code{\link{xcountCutoff}} to determine whether the number of tables to examine is greater than \code{cutoff}. If it is, then \code{\link{mtest}} is used. Otherwise \code{\link{xtest}} is used. The result is a robust test which will always provide a meaningful and accurate P value. Each table examined is compared with the observed counts according to each of four measures of fit:  \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} corresponding to the log-likelihood ratio, the null-hypothesis probability, the U-score or the Pearson X^2 value. It can also plot a histogram showing the distribution of any of these statistics.
#' 
#' 
#' @references The methods are described by \href{http://dx.doi.org/10.1534/genetics.109.108977}{Engels, 2009. \bold{Genetics} 183:1431}.
#' 
#' @importFrom parallel mclapply
#' 
#' 
#' @param c The genotype counts. You must provide the number of each genotype. So if there are \eqn{k} alleles, you need to include the number of each of the \eqn{k(k+1)/2} genotypes. The format of \code{x} is somewhat flexible: It can be a square matrix, but only the lower-left half is used. It can be a vector of the observations in the order \eqn{a_11, a_21, a_22, a_31, ..., a_kk}. For compatability with the packages \code{genetics} and \code{adegenet}, it can also be an object of class \code{genind}, \code{genotype}, or a \code{data.frame}. If \code{c} contains multiple samples, the \code{parallel} package will be used in an attempt to employ multi-cores.
#' @param method Can be \dQuote{auto}, \dQuote{exact} or \dQuote{monte} to indicate the method to use. If \dQuote{auto}, the \code{hwx.test} will first check to see whether the total number of tables exceeds a cutoff specified by the parameter \code{cutoff}.
#' @param cutoff If \code{method} is set to \dQuote{auto}, then \code{cutoff} is used to decide whether to perform the test via the full enumeration or Monte Carlo method. If the number of tables is less than \code{cutoff}, then a full enumeration is performed. Otherwise the method will be Monte Carlo with \code{B} random trials.
#' @param B The number of trials to perform if Monte Carlo method is used
#' @param statName can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.
#' @param histobins If 0, no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If \code{histobins} is set to a number greater than 1, a histogram with \code{histobins} bins is plotted.
#' @param histobounds A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{hwx.test} will compute reasonable bounds to include most of the distribution.
#' @param showCurve whether to show a blue curve indicating the asymptotic (chi squared) distribution. This only works for \code{LLR} and \code{Chisq}
#' @param safeSecs After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.
#' @param detail Determines how much detail is printed. If it is set to 0, nothing is printed (useful if you use \code{hwx.test} programmatically.)
#' 

#' 
#' @return Returns a list of class \code{hwtest} which includes the following items:
#' \item{$ Pvalues}{The four computed P values corresponding to the test statistics: \code{LLR}, \code{Prob}, \code{U} and \code{Chisq} in that order.}
#' \item{$ observed}{The four observed statistics in the same order as above}
#' \item{$ ntrials}{The number of tables examined during the calculation if done by Monte Carlo}
#' \item{$ tableCount}{The total number of tables if done by full enumeration}
#' \item{$ genotypes}{The input matrix of genotype counts}
#' \item{$ alleles}{The allele counts \eqn{m} corresponding to the input genotype counts}
#' \item{$ statName}{Which statistic to use for the histogram and in the \code{p.value} item}
#' \item{$ method}{Which method was used, \dQuote{exact} or \dQuote{monte}}
#' \item{$ detail}{An integer indicating how much detail to print. Use 0 for no printing}
#' \item{$ SE}{vector with the standard error for each stat. Only applicable with Monte Carlo tests}

#' 
#' @examples
#' # Data from Louis and Dempster 1987 Table 2 and Guo and Thompson 1992 Figure 2:
#' c <- c(0,3,1,5,18,1,3,7,5,2)
#' hwx.test(c)
#' # To see a histogram of the LLR statistic:
#' hwx.test(c, histobins=TRUE)
#' # For a histogram of the U statistic and other details of the result:
#' hwx.test(c, statName="U", histobins=TRUE, detail=3)
#' 

#' @export
hwx.test <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	UseMethod("hwx.test")
}



#' @export
hwx.test.matrix <- 
function(c,  method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	nAlleles <- length(alleleCounts(c));
	if(nAlleles < 2) return(hwx.test.logical(FALSE))
	statNames <- c("LLR", "Prob", "U", "Chisq");
	statID <- which(statNames==statName);
	c <- remove.missing.alleles(c)
	if(method=="auto") method <- if(xcountCutoff(c, cutoff))method <- "monte" else method <- "exact"
	if(method=="monte") 
		value <- mtest(c, ntrials=B, statName, histobins, histobounds, showCurve, safeSecs, detail)
	else
		value <- xtest(c, statName, histobins, histobounds, showCurve, safeSecs, detail);
	value$method=method;
#	value$p.value=value$Pvalues[statID];
	value$statName=statName;
	value$detail=detail;
	class(value) <- "hwtest"
	return(value)
}

#' @export
hwx.test.integer <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	c <- vec.to.matrix(c);
	hwx.test.matrix(c,  method=method, cutoff=cutoff, B=B, statName=statName, histobins=histobins, histobounds=histobounds, showCurve=showCurve, safeSecs=safeSecs, detail=detail)
}
#' @export
hwx.test.numeric <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	c <- vec.to.matrix(c);
	hwx.test(c,  method=method, cutoff=cutoff, B=B, statName=statName, histobins=histobins, histobounds=histobounds, showCurve=showCurve, safeSecs=safeSecs, detail=detail)
}

#' @export
hwx.test.table <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	hwx.test(unclass(c),method=method, cutoff=cutoff, B=B, statName=statName, histobins=histobins, histobounds=histobounds, showCurve=showCurve, safeSecs=safeSecs, detail=detail)
}

#' @export
hwx.test.genotype <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	requireNamespace("genetics")
	tab <- table(factor(genetics::allele(c, 1), levels=genetics::allele.names(c)), factor(genetics::allele(c, 2), levels=genetics::allele.names(c)));
	hwx.test(unclass(t(tab)), method=method, cutoff=cutoff, B=B, statName=statName, histobins=histobins, histobounds=histobounds, showCurve=showCurve, safeSecs=safeSecs, detail=detail)
}


#' @export
hwx.test.genind <- function(c, method = "auto", cutoff = 1e+07, B = 1e+05, statName = "LLR", histobins = 0, histobounds = c(0, 0), showCurve = T, safeSecs = 100, detail = 2) {
	if (requireNamespace("adegenet")) {
		if (!(adegenet::is.genind(c))) 
			stop("function requires a genind object")
		if (c@ploidy[1] != as.integer(2)) 
			stop("function requires diploid data")
		df <- adegenet::genind2df(c, pop = c@pop, sep = "/")
		outcome <- hwx.test(df, method = method, cutoff = cutoff, B = B, statName = statName, histobins = histobins, histobounds = histobounds, showCurve = showCurve, safeSecs = safeSecs, detail = detail)
	} else {
		cat("\nERROR: package adegenet is not installed.\n")
		outcome <- hwx.test.logical(c)
	}
	return(outcome)
}


#' @export
hwx.test.data.frame <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	mlist <- df.to.matrices(c);
	hwx.test(mlist,method=method, cutoff=cutoff, B=B, statName=statName, histobins=histobins, histobounds=histobounds, showCurve=showCurve, safeSecs=safeSecs, detail=detail)
}

#' @export
hwx.test.list <- function(c, method = "auto", cutoff = 1e+07, B = 1e+05, statName = "LLR", histobins = 0, histobounds = c(0, 0), showCurve = T, safeSecs = 100, 
	detail = 2) {
	nameItems <- function(x, nam = "") {
		if (class(x) == "hwtest") 
			x$sampleName <- nam
		if (class(x) == "list") {
			for (na in names(x)) x[[na]] <- nameItems(x[[na]], paste(na, nam))
		}
		return(x)
	}

	cores <- getOption("mc.cores", 1L)
	if (cores >= 1 && requireNamespace("parallel") && Sys.info()[1] != "Windows") {
		RNGkind("L'Ecuyer-CMRG")
		resultList <- parallel::mclapply(c, hwx.test, method = method, cutoff = cutoff, B = B, statName = statName, histobins = histobins, histobounds = histobounds, 
			showCurve = showCurve, safeSecs = safeSecs, detail = detail, mc.allow.recursive = T, mc.cores = cores)
	} else {
		resultList <- lapply(c, hwx.test, method = method, cutoff = cutoff, B = B, statName = statName, histobins = histobins, histobounds = histobounds, 
			showCurve = showCurve, safeSecs = safeSecs, detail = detail)
	}
	return(nameItems(resultList))
}


#' @export
hwx.test.logical <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	statNames <- c("LLR", "Prob", "U", "Chisq");
	a <- list(Pvalues=rep(NA, 4), observed=rep(NA, 4), method=NA, statName=NA);
	names(a$Pvalues) <- statNames;
	names(a$observed) <- statNames;
	class(a) <- "hwtest"
	return(a)
}