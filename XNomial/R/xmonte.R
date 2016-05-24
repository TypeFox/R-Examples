# xmonte
# Goodness-of-fit test for multinomial distribution with fixed probabilities
# Monte Carlo version
# Bill Engels
#
# Compute the probability that the test statistic is "at least as extreme" as observed. The test statistic can
# be selected from among three choices: "LLR" for Log Likelihood Ratio, "Prob" for the probability itself or "Chisq"
# for the standard chi square statistic.

#' Perform Multinomial Goodness-Of-Fit Test By Monte-Carlo Simulations
#' 
#' Use \code{xmonte} to compute a P value to test whether a set of counts fits a specific multinomial distribution. It does this by examining a large number of random outcomes and finding the probability of those cases which deviate from the expectation by at least as much as the observed.
#' 
#' @inheritParams xmulti
#' @param ntrials the number of random trials to look at, such as \code{ntrials=100000}
#' 
#' @return \code{xmonte} returns a list with the following components:
#' \item{$ obs}{the observed numbers used as imput}
#' \item{$ expr}{expected ratios, arbitrary scale}
#' \item{$ ntrials}{the number of random tables examined}
#' \item{$ statType}{which test statistic was used}
#' \item{$ pLLR/pProb/pChi}{the P value computed for the given test statistic}
#' \item{$ standard.error}{the binomial standard error of the estimated P value}
#' \item{$ observedLLR}{the value of LLR statistic for these data}
#' \item{$ observedProb}{the multinomial probability of the observed data under the null hypothesis}
#' \item{$ observedChi}{observed value of the chi square statistic}
#' \item{$ histobins}{number of bins in the histogram (suppressed if zero)}
#' \item{$ histobounds}{range in histogram (suppressed if not used)}
#' \item{$ histoData}{data for histogram (suppressed if not used) Length is \code{histobins}}
#' \item{$ asymptotoc.p.value}{the P value obtained from the classical asymptotic test -- use for comparison only}
#' 
#' @examples
#' #One of Gregor Mendel's crosses produced four types of pea seeds in the numbers as follows:
#' peas <- c(315, 108, 101, 32)
#' #and he expected them to appear in the ratio of 9:3:3:1 according to his genetic model.
#' expected <- c(9, 3, 3, 1)
#' #Test Mendels theory using
#' xmonte(peas, expected)
#' #To see a histogram of the likelihood ratio statistic, use:
#' xmonte(peas, expected, histobins = TRUE)
#' #The red areas of the histogram represent those outcomes deviating from the expected 9:3:3:1 ratio 
#' #at least as much as the observed numbers. (Much has been made of the tendency for Mendel's data 
#' #to fit the expectations better than expected!)
#' #If you wish to use the standard chisquare statistic as a measure of goodness-of-fit instead 
#' #of the LLR, use:
#' xmonte(peas, expected, statName="Chisq", histobins=TRUE)

#' @useDynLib XNomial
#' @export
xmonte <-
function(obs, expr, ntrials = 100000, statName = "LLR", histobins = F, histobounds = c(0,0), showCurve = T, detail = 1, safety=1e8) {
  		if(length(obs) != length(expr)) stop("\nThere must be the same number of observed and expected values\n");
		if(length(obs) < 2) stop("\nThere must be at least two categories\n");
  		if(ntrials > safety) stop("\nTo run this many trials, you must increase the 'safety' parameter\n")
  		statType = which(c("LLR", "Prob", "Chisq") == statName);
  		if(length(statType) == 0) {statType <- 4}
  		if(histobins == 1){ histobins  <- 500};   #The default is 500 bins
  		if(histobounds[[2]] == 0) {histobounds[[2]] <- qchisq(.999, length(obs) - 1)}
  		invisible(rmultinom(4, obs, expr)); # just to initialize the rmultinom function
   		result = .C("montenomialTest",
		  	obs=as.integer(obs),
		  	expr=as.double(expr),
		  	ntrials = as.integer(ntrials),
		  	nn = as.integer(length(obs)),
		  	statType=as.integer(statType), #1 for LLR, 2 for prob, 3 for chi
		  	pLLR = as.double(0),
		  	pProb = as.double(0),
		  	pChi = as.double(0),
		  	observedLLR = as.double(logLRmultinomial(obs, expr)),
		  	observedProb = as.double(dmultinom(obs, prob = expr)),
		  	observedChi = as.double(chiStat(obs, expr)),
		  	histobins = as.integer(histobins),
		  	histobounds = as.double(histobounds),
		  	histoData = as.integer(integer(if(histobins) {histobins} else {1}))
		  	,PACKAGE = "XNomial"
		  	);
		pval <- switch(statType, result$pLLR, result$pProb, result$pChi);
		standardError <- sqrt(pval*(1-pval)/ntrials); #binomial standard error
		if(detail >= 1) {cat("\nP value (",statName, ") = ", pval, " \u00b1 ", formatC(standardError), sep="")};
		# \u00b1 is the plus-or-minus symbol
		if(detail >= 2) {cat("\n", formatC(ntrials), " random trials\n Observed: ", obs, "\n Expected Ratio: ", expr, sep=" ")};
		if(histobins){statHistogram(result, showCurve)}
		value  <- result;
		value[["nn"]] <- NULL;
	    if(histobins==0){   #suppress histogram information if not used
		  value[["histobins"]] <- NULL;
		  value[["histobounds"]] <- NULL;
		  value[["histoData"]] <- NULL;
		  }
		#include only the relevant P value
		if(statType != 1) {value[["pLLR"]] <- NULL};
		if(statType != 2) {value[["pProb"]] <- NULL};
		if(statType != 3) {value[["pChi"]] <- NULL}
		#add some other info that might be useful
	    value$statType <- statName;
	    asymptotic.p.value  <-  pchisq(result$observedChi, result$nn-1, lower.tail=F);
	    value <- c(value, asymptotic.p.value = asymptotic.p.value);
   	    value <- append(value, c(standard.error=standardError), after=which(names(value)=="statType")+1)
	    invisible(value);
  	}
