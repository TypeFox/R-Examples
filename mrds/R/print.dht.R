#' Prints density and abundance estimates based on Horvitz-Thompson-like estimator
#'
#' Outputs summary statistics, abundance and density by region (if any) and optionally a correlation matrix if more than one region.
#'
#' @method print dht
#' @param x dht object that results from call to dht for a specific ddf object
#' @param cor if TRUE outputs correlation matrix of estimates
#' @param bysample if TRUE, prints results for each sample
#' @param vcmatrices if TRUE, prints variance-covariance matrices
#' @param \dots unspecified and unused arguments for S3 consistency
#' @export
#' @return None
#' @author Jeff Laake
#' @seealso \code{\link{dht}}
#' @keywords utility
print.dht <- function(x,cor=FALSE,bysample=FALSE,vcmatrices=FALSE,...){

  print.tables <- function(x,cor,bysample,vcmatrices){
    cat("\nSummary statistics:\n")
    print(x$summary)
    if("N" %in% names(x)){
      cat("\nAbundance:\n")
      print(x$N)
    }
    cat("\nDensity:\n")
    print(x$D)
    if(vcmatrices){
      cat("\nAbundance variance-covariance matrix component from estimating detection function\n" )
      print(x$vc$detection$variance)
      cat("\nAbundance variance-covariance matrix component from sample selection\n" )
      print(x$vc$er)
    }
    if(cor){
      cat("\nCorrelation matrix:\n")
      print(x$cormat)
    }
    if(bysample){
      cat("\nEstimates by sample:\n")
      print(x$bysample)
    }
  }

  if(is.null(x$clusters)){
    print.tables(x$individuals,cor,bysample,vcmatrices)
  }else{
    cat("\nSummary for clusters\n")
    print.tables(x$clusters,cor,bysample,vcmatrices)
    cat("\nSummary for individuals\n")
    print.tables(x$individuals,cor,bysample,vcmatrices)
    cat("\nExpected cluster size\n")
    #Added CV as an output LJT 14/10/09
    S <- x$Expected.S
    if(!is.null(S$se.Expected.S)){
      S$cv.Expected.S <- S$se.Expected.S/S$Expected.S
      S$cv.Expected.S[S$Expected.S==0] <- 0
    }
    print(S)
  }
  invisible()
}
