#' Random generation of missingness in a data frame.
#'
#' Generates random \code{NA}'s in in a column or groups of columns of a data frame.
#' Used in imputation simulations based on complete datasets.        
#'
#' @param aDataFrame A data.frame where missingness will be simulated.
#' @param columnNames A string or a vector of strings that describe the column 
#'    names (labels) where missingness will be simulated.
#' @param percentMissing The percentage of missingness within specified columns.
#'    "Percent missing" uses a binomial distribution to simulate missing data.
#'    Default is 10 (i.e. 10\% missing).  Use \code{\link{impute_missingness}} for
#'    a summary of these randomly generated missing data.
#'
#' @return A data table with columns of missing data (specified as \code{NA}'s).
#'
#' @importFrom stats rbinom
#' @export random_missingness

random_missingness <- function(aDataFrame, 
                               columnNames, 
                               percentMissing = 10) {
                               
  missingness <- !rbinom(nrow(aDataFrame), 1, 1.0 - (percentMissing/100))
  for(i in columnNames) {
    newColumn <- as.matrix(aDataFrame[i])
    newColumn[missingness] <- NA
    aDataFrame[i] <- newColumn
  }
  return(aDataFrame)
}
