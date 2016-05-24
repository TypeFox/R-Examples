#' Convert output of ibs.pairwise.db for plotting using DNAtools package
#' 
#' @param m Matrix with the number of full/partial matches (output of \code{\link{ibs.pairwise.db}})
#' @details Converts the matrix to a \code{dbcompare} object for use with the plot function of the \code{DNAtools} package.
#' @return Object of class \code{dbcompare} (compatible with \code{DNAtools} package)
#' @seealso \link{ibs.pairwise.db}
#' @examples
#' data(freqsNLsgmplus)
#' db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
#' M <- as.dbcompare(ibs.pairwise.db(db))
#' \dontrun{
#' require(DNAtools)
#' plot(M)
#' }
#' @export
as.dbcompare <- function(m){
  ret <- list(m=m$M,excludedProfiles="none")
  attributes(ret)$call <- list(loci=(ncol(m$M)-1),single=0,collapse=FALSE,vector=FALSE,wildcard=c(FALSE,FALSE))
  class(ret) <- "dbcompare"
  
  ret
}
NULL