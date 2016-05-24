##
## Tests based on adjusted residual
##
## --------------------------------------------------------------------------

rvgt.Mtest <- function(ftable)

  ## ------------------------------------------------------------------------
  ## Perform M-test on frequency table.
  ## ------------------------------------------------------------------------
  ## ftable : Object of class "rvgt.ftable" containing frequencies
  ## ------------------------------------------------------------------------
  ## Return:
  ## list of p values of cumulative frequencies
  ## ------------------------------------------------------------------------
  ## References:
  ## [1] C. Fuchs and R. Kenett. A test for Detecting Outlying Cells in
  ## the Multinomial Distribution and Two-Way Contingency Tables
  ## Journal of American Statistical Association, Vol 75, Jun 1980, 395-398
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
  
  ## vector to store p-values.
  pval <- numeric(r)
  
  ## array for computing cumulative frequencies
  fcum <- numeric(nbins)
  
  ## array for storing calculated M statistics
  maxz <- numeric(r)

  ## Iteration for calculating p-values
  for(i in 1:r){
      fcum <- fcum + table[i,]
      
      ## total samplesize
      nt <- i * n

      ## calculating adjusted residual
      z <- (fcum - nt*p0) / sqrt(nt*p0*(1-p0))
           
      ## Calculating absolute maximum of z'j
      maxz<-max(abs(z))

      ## Calculating p-value
      pval[i] <- min(1, 2*nbins*(1-pnorm(maxz)))
      pval[i] <- max(pval[i],1.e-300)
  }   
  
  ## return result as object of class "rvgt.htest"
  result <- list (type="M-test", n=n, rep=r, breaks=nbins+1, pval=pval)
  class(result) <- "rvgt.htest"

  return (result)    
} 

##-------------------------------------------------------------------------
