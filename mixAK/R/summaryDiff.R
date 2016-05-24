##
##  PURPOSE:  Function to calculate (posterior) summary
##            statistics for a difference of two quantities supplied as (MCMC) samples
##            * primarily used to calculate posterior summary for the difference of the
##              deviances of two competing models 
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   19/05/2010 (as a stand alone function)
##             26/11/2010:  added to the mixAK package
##
##  FUNCTIONS: summaryDiff
##
## ==========================================================================

## *************************************************************
## summaryDiff
## *************************************************************
##
summaryDiff <- function(x, y, prob=c(0.025, 0.5, 0.975), cut=c(-2*log(9), 0), na.rm=TRUE)
{
  if (any(prob < 0)) stop("all prob values must be nonnegative")
  if (any(prob > 1)) stop("all prob values must not exceed 1")  

  #if (missing(n)){
    x <- as.numeric(x)
    y <- as.numeric(y)
  #}else{
  #  x <- as.numeric(x)[sample.int(length(x), size=n, replace=FALSE)]
  #  y <- as.numeric(y)[sample.int(length(y), size=n, replace=FALSE)]    
  #}  

  if (length(x) != length(y)) stop("x and y must be of the same length")
  diff <- x - y
  
  #diff <- rep(x, each=length(y)) - rep(y, length(x))                    ### all pairwise differences

  Ediff <- mean(diff[!(diff == Inf | diff == -Inf)], na.rm=na.rm)
  #Ediff <- mean(diff, na.rm=na.rm)
  Qdiff <- quantile(diff, prob=prob, na.rm=na.rm)
  Pcut <- numeric(length(cut))
  for (i in 1:length(cut)){
    ptab <- prop.table(table(diff < cut[i]))
    if ("TRUE" %in% names(ptab)) Pcut[i] <- ptab["TRUE"] else Pcut[i] <- 0
  }  
  names(Pcut) <- paste("P(diff < ", round(cut, 2), ")", sep="")
  summ <- c(Ediff, Qdiff)
  names(summ)[1] <- "Mean"
  
  RET <- list(summary=summ, Pcut=Pcut)
  return(RET)    
}  
