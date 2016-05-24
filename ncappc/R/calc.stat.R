# R function to calculate a set of statistics for a given array of numbers
# Chayan, 12/2014

# roxygen comments
#' Calculates a set of statistics for a given array of numbers.
#'
#' \pkg{calc.stat} calculates a set of statistics for a given array of numbers.
#'
#' \pkg{calc.stat} calculates a set of statistics for a given array of numbers.
#' The calculated statistics are
#' \itemize{
#'  \item Ntot = length of the array
#'  \item Nunique = Number of unique elements
#'  \item Min = Minimum value of the array
#'  \item Max = Maximum value of the array
#'  \item Mean = Mean value of the array
#'  \item Median = Median value of the array
#'  \item SD = Standard deviation value of the array
#'  \item SE = Standard error value of the array
#'  \item CVp = Percent coefficient of variation of the array
#'  \item CI95u = Upper limit of the 95\% confidence interval of the array
#'  \item cI95l = Lower limit of the 95\% confidence interval of the array
#'  \item gMean = Geometric mean value of the array
#'  \item gCVp = Geometric percent coefficient of variation of the array
#' }
#' 
#' @param x a numeric array
#'
#' @return An array of calculated statistics of a given set of numbers
#' @export
#'

calc.stat <- function(x){
  "median" <- "sd" <- "qt" <- "tail" <- "head" <- "lm" <- "coef" <- NULL
  rm(list=c("median","sd","qt","tail","head","lm","coef"))
  x       <- as.numeric(x)
  Ntot    <- length(x)
  Nunique <- length(unique(x))
  Min     <- min(x)
  Max     <- max(x)
  Mean    <- mean(x)
  Median  <- median(x)
  SD      <- sd(x)
  SE      <- SD/sqrt(Ntot)
  CVp     <- 100*SD/Mean
  CI95l   <- Mean-abs(qt(0.025,Ntot-1)*SE)  # two-tailed
  CI95u   <- Mean+abs(qt(0.975,Ntot-1)*SE)  # two-tailed
  
  gm_mean = function(x, na.rm=TRUE){
    if(any(x < 0, na.rm = TRUE)){
      return(NaN)
    }
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  }
  
  gm_cv_p = function(x, na.rm=TRUE){
    if(any(x <= 0, na.rm = TRUE)){
      return(NaN)
    }
    sqrt(exp(sd(log(x))^2)-1)*100
  }
  
  gMean   <- gm_mean(x)
  gCVp    <- gm_cv_p(x)
  stPrm   <- c(Ntot,Nunique,Min,Max,Mean,Median,SD,SE,CVp,CI95l,CI95u,gMean,gCVp)
  names(stPrm) <- c("Ntot","Nunique","Min","Max","Mean","Median","SD","SE","CVp","CI95l","CI95u","gMean","gCVp")
  return(stPrm)
}


