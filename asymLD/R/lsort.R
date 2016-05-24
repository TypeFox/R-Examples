#' Sort a data.frame.
#'
#' A function to sort a data.frame on specific columns. 
#' 
#' @param dat	a dataframe or a matrix ("dimnames" are used as the variable names for a matrix)
#' @param by	a vector or a list of variable names or column indices specifying the "sort by" 
#'   variables, the default is to sort by all variables in the order they appear in the data set.
#' @param asc	a vector with the same length as "by" indicating whether the sorting of each "by" 
#'   variable is in ascending order, the default order is ascending.
#' @param na.last	a flag indicating whether missing values are placed as the last 
#'   elements in the data set, the default is TRUE
#'
#' @return The return value is a sorted dataframe. 
#'
#' @section Details:
#' The input dataframe is not modified. The code is adapted from code posted to an old s-news 
#' listserve.
#'
#' @examples
#' \dontrun{
#' library(asymLD)
#' data(snp.freqs)
#' 
#' # sort snp.freqs by "locus1" (ascending) and "allele1" (descending)
#' newdata <- lsort(snp.freqs, by=c("locus1","allele1"), asc=c(T,F))
#' head(newdata)
#' # sort snp.freqs by the fourth and the second variable (ascending)
#' newdata <- lsort(snp.freqs, by=c(4,2))
#' # sort "snp.freqs" by "locus1" and the 5th variable (ascending)
#' newdata <- lsort(snp.freqs, by=list("locus1",5))
#' 
#' }
#' @export

lsort <- function(dat, by = 1:dim(dat)[2], asc = rep(TRUE, length(by)), na.last = TRUE) {
  m <- dim(dat)[1]
  keys <- 1:m
  rotate <- m:1
  for(i in length(by):1) {
          if(asc[i])
               keys[] <- keys[sort.list(dat[, by[[i]]][keys], na.last= na.last)]
          else keys[] <- keys[order(dat[, by[[i]]][keys], rotate, na.last = na.last)[rotate]]
  }
  dat[keys,  ]
}


