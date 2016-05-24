##
## Tests based on chisquare goodness-of-fit tests
##
## --------------------------------------------------------------------------

rvgt.chisq <- function (ftable)

  ## ------------------------------------------------------------------------
  ## Perform chisquare test on frequency table.
  ## ------------------------------------------------------------------------
  ## ftable : Object of class "rvgt.ftable" containing frequencies
  ## ------------------------------------------------------------------------
  ## Return:
  ## list of p values of cumulative frequencies
  ## ------------------------------------------------------------------------
{ 
  ## check arguments
  if (missing(ftable) || class(ftable) != "rvgt.ftable")
    stop ("Argument 'ftable' missing or invalid.")

  ## get table
  table <- ftable$count
  
  ## samplesize for one repetition
  n <- ftable$n

  ## number of repetitions
  r <- ftable$rep
  
  ## number of bins
  nbins <- ncol(table)

  ## probabilities under null hypothesis
  ubreaks <- ftable$ubreaks
  p0 <- diff(ubreaks)
  
  ## Vector to store p-values.
  pval <- numeric(r)
  
  ## array for computing cumulative frequencies
  fcum <- numeric(nbins)

  ## compute p-values of cumulative frequencies
  for (i in 1:r) {
    fcum <- fcum + table[i,]
    pval[i] <- chisq.test(fcum,p=p0)$p.value
  }
  
  ## return result as object of class "rvgt.htest"
  result <- list (type="chisq", n=n, rep=r, breaks=nbins+1, pval=pval)
  class(result) <- "rvgt.htest"

  return (result)    
}

## --------------------------------------------------------------------------

## EXPERIMENTAL code
## Not exported!

rvgt.chisq.level2 <- function (ftable)

  ## ------------------------------------------------------------------------
  ## Perform a level-2 chisquare test on frequency table.
  ## level-1: chisquare test on each row of table
  ## level-2: test p-values for uniformity using test by Fisher (1967)
  ## ------------------------------------------------------------------------
  ## ftable : Object of class "rvgt.ftable" containing frequencies
  ## ------------------------------------------------------------------------
  ## Return:
  ## list of p-values of cumulative frequencies
  ## ------------------------------------------------------------------------
  ## References:
  ## [1] R. A. Fisher (1967). Statistical Methods for Research Workers.
  ##     4th edition, New York.
  ## [2] M. A. Stephens (1986). Tests for the uniform distribution, p.358.
  ##     in: R. B. D'Agostino and M. A. Stephens (eds.), Goodness-of-fit
  ##     techniques, New York: Dekker.
  ## ------------------------------------------------------------------------
  ## REMARK: This function is EXPERIMENTAL and yet not exported!
  ## ------------------------------------------------------------------------
{ 
  ## check arguments
  if (missing(ftable) || class(ftable) != "rvgt.ftable")
    stop ("Argument 'ftable' missing or invalid.")

  ## get table
  table <- ftable$count
  
  ## samplesize for one repetition
  n <- ftable$n

  ## number of repetitions
  r <- ftable$rep
  
  ## number of bins
  nbins <- ncol(table)
 
  ## Vector to store p-values.
  pval <- numeric(r)
  ptmp <- numeric(r)

  ## compute level-2 p-values for increasing number of repetitions
  for (i in 1:r) {
    ## p-value for i-th sample
    ptmp[i] <- chisq.test(table[i,])$p.value
    ## test statistic (chisq distributed with 2*i degrees of freedom)
    P2 <- -2 * sum(log(1-ptmp[1:i]))
    ## pval = Prob(X^2 <= P2) 
    pval[i] <- pchisq( q=P2, df=2*i, lower.tail = TRUE)
  }

  ## return result as object of class "rvgt.htest"
  result <- list (type="chisq-2", n=n, rep=r, breaks=nbins+1, pval=pval)
  class(result) <- "rvgt.htest"

  return (result)
}

## --------------------------------------------------------------------------
