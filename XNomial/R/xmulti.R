# xmulti
# Goodness-Of-Fit Test For Multinomial Distribution With Fixed Probabilities
# Full Enumeration version
# Bill Engels
#Updated to trap underflow errors 16 December 2015
#
# Compute the probability that the test statistic is "at least as extreme" as observed. The test statistic can
# be selected from among three choices: \dQuote{LLR} for Log Likelihood Ratio, \dQuote{Prob} for the probability itself or \dQuote{Chisq}
# for the standard chi square statistic. 

# note: n * log(pmax) must not be less than .Machine$double.xmin or there will be underflow.

#' Perform Multinomial Goodness-Of-Fit Test By Full Enumeration
#' 
#' Use \code{xmulti} to compute a P value to test whether a set of counts fits a specific multinomial distribution. It does this by examining all possible outcomes with the same total count and determining the total (multinomial) probability of those cases which deviate from the expectation by at least as much as the observed. Please see the vignette for more.
#' 
# These calls to eliminate the "no visible global function definition for ..." from CRAN submission
#' @importFrom graphics barplot lines
#' @importFrom stats dchisq dmultinom pchisq qchisq rmultinom
#' 

#' @param obs vector containing the observed numbers. All are non-negative integers summing to \code{> 0}.
#' @param expr vector containing expectation. The length should be the same as that of \code{obs} and they should be non-negative summing to \code{> 0}. They need not be integers or sum to one.
#' @param statName name of the test statistic to use as a measure of how deviant an observation is from the expectation. The choices are: \dQuote{LLR} for the log-likelihood ratio, \dQuote{Prob} for the probability, \dQuote{Chisq} for the chisquare statistic.
#' @param histobins specifies histogram plot. If set to 0, \code{F} or \code{FALSE} no histogram is plotted. If set to 1 or \code{T} or \code{TRUE} a histogram with 500 bins will be plotted. If set to a number \code{> 1} a histogram with that number of bins is plotted.
#' @param histobounds vector of length 2 indicating the bounds for the histogram, if any. If unspecified, bounds will be determined to include about 99.9 percent of the distribution.
#' @param showCurve should an asymptotic curve be drawn over the histogram?
#' @param detail how much detail should be reported concerning the P value. If 0, nothing is printed for cases where the function is used programmatically. Minimal information is printed if \code{detail} is set to 1, and additional information if it is set to 2.
#' @param safety a large number, such as one billion, to set a limit on how many samples will be examined. This limit is there to avoid long computations.
#' 
#' @return \code{xmulti} returns a list with the following components:
#' \item{$ obs}{the observed numbers used as imput}
#' \item{$ expr}{expected ratios}
#' \item{$ statType}{which test statistic was used}
#' \item{$ pLLR}{the P value with LLR as the test statistic}
#' \item{$ pProb}{the P value with the multinomial probability as test statistic}
#' \item{$ pChi}{the P value with the chisquare as test statistic}
#' \item{$ observedLLR}{the value of LLR statistic for these data}
#' \item{$ observedProb}{the multinomial probability of the observed data under the null hypothesis}
#' \item{$ observedChi}{observed value of the chi square statistic}
#' \item{$ histobins}{number of bins in the histogram (suppressed if zero)}
#' \item{$ histobounds}{range in histogram (suppressed if not used)}
#' \item{$ histoData}{data for histogram (suppressed if not used) Length is \code{histobins}}
#' \item{$ asymptotoc.p.value}{the P value obtained from the classical asymptotic test}
#' \item{$ cases.examined}{the total number of possible tables}

#' 
#' @examples
#' #
#' #One of Gregor Mendel's crosses produced four types of pea seeds in the numbers:
#' #
#' peas <- c(315, 108, 101, 32)
#' #
#' #and he expected them to appear in the ratio of 9:3:3:1 according to his genetic model.
#' #
#' expected <- c(9, 3, 3, 1)
#' #
#' #Test Mendel's theory using
#' #
#' xmulti(peas, expected)
#' #
#' #In this example, the number of cases examined was 28956759,
#' #and it probably took your computer less than half a second.
#' #To see a histogram of the likelihood ratio statistic, use:
#' #
#' xmulti(peas, expected, histobins = TRUE)
#' #
#' #The red areas of the histogram represent those outcomes deviating from the expected 9:3:3:1 ratio 
#' #at least as much as the observed numbers. (Much has been made of the tendency for Mendel's data 
#' #to fit the expectations better than expected!)
#' #If you wish to use the standard chisquare statistic as a measure of goodness-of-fit instead 
#' #of the LLR, use:
#' #
#' xmulti(peas, expected, statName="Chisq", histobins=TRUE)


#' @useDynLib XNomial
#' @export
xmulti <-
function(obs, expr, statName = "LLR", histobins = F, histobounds = c(0,0), showCurve = T, detail=1, safety=1e9) {
		if(length(obs) != length(expr)) stop("\nThere must be the same number of observed and expected values\n");
		if(length(obs) < 2) stop("\nThere must be at least two categories\n");
		ntotal <- sum(obs)
		extotal <- sum(expr)
		if(ntotal <= 0) stop("\nThe total observations must be positive")
		if(extotal <= 0) stop("\nThe sum of the expected numbers must be positive")
		maxprob  <- max(expr/extotal)
		if (ntotal * log(maxprob) <= log(.Machine$double.xmin)) stop("\nThe numbers of observations are too large causing underflow error. The monte carlo version, \"xmonte\" is recommended for this case.")
  		statType = which(c("LLR", "Prob", "Chisq") == statName);
  		if(histobins == 1){ histobins  <- 500};   #The default is 500 bins
  		if(histobounds[[2]] == 0) {histobounds[[2]] <- qchisq(.999, length(obs) - 1)};
  		ntrials  <- ntables(obs);
  		if(ntrials > safety) {
  			cat("Full enumeration requires examination of", ntrials, "tables.\n")
  			# \u2022 is a bullet symbol
  			stop("This operation could take more than several minutes.\n    \u2022 The monte carlo version, \"xmonte\" is recommended for this case. \n    \u2022 To override this cutoff, change the parameter 'safety' to something greater than the required number of trials.")
  		}
  		result = .C("exactMultinomialTest",
		  	obs=as.integer(obs),
		  	expr=as.double(expr),
		  	nn=as.integer(length(obs)),
		  	statType=as.integer(statType), #0 for LLR, 1 for prob
		  	pLLR = as.double(0),
		  	pProb = as.double(0),
		  	pChi = as.double(0),
		  	observedLLR = as.double(logLRmultinomial(obs, expr)),
		  	observedProb = as.double(dmultinom(obs, prob = expr)),
		  	observedChi = as.double(chiStat(obs, expr)),
		  	histobins = as.integer(histobins),
		  	histobounds = as.double(histobounds),
		  	histoData = as.double(double(if(histobins) {histobins} else {1}))
		  	,PACKAGE = "XNomial"
		  	);
		pval  <- c(result$pLLR, result$pProb, result$pChi)[statType];  	
	  	if(detail == 1){
			cat("\nP value (", statName, ") = ", formatC(pval), "\n", sep=""); 		
	  	}
	  	if(detail >= 2) {
	  		cat("\nP value  (LLR)  = ", formatC(result$pLLR));
	  		cat("\nP value (Prob)  = ", formatC(result$pProb));
	  		cat("\nP value (Chisq) = ", formatC(result$pChi));
	  	}
	  	if(detail >= 3) {
	  		cat("\n\nObserved: ", obs, "\nExpected ratio: ", expr, "\nTotal number of tables: ", ntrials,"\n")
	  	}
		  if(histobins){statHistogram(result, showCurve)}
	  value  <- result;
	  value[["nn"]] <- NULL;
	  if(histobins==0){
		  value[["histobins"]] <- NULL;
		  value[["histobounds"]] <- NULL;
		  value[["histoData"]] <- NULL;
	  }
	  value$statType <- statName;
	  asymptotic.p.value  <-  pchisq(result$observedChi, result$nn-1, lower.tail=F);
	  value <- c(value, asymptotic.p.value = asymptotic.p.value, cases.examined=ntrials);
	  invisible(value);
  	}
