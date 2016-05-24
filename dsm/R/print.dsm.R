#' Print a description of a density surface model object
#'
#' This method just gives a short description of the fitted model. Use the
#' \code{\link{summary.dsm}} method for more information.
#'
#' @param x a \code{dsm} object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @export
#' @author David L. Miller
#' @seealso \code{\link{summary.ds}}
#' @keywords utility
print.dsm<-function(x,...){

  # the code here is chopped together from mgcv and mrds

  # if we fit a gamm then we should just grab the gam bit for this
  # but note that a gamm was used!
  gamm <- FALSE
  if("gamm" %in% class(x)){
    x<-x$gam
    gamm <- TRUE
  }

  ### General information

  cat("\nDensity surface model\n")
  cat("Response : ",as.character(x$formula)[2] , "\n")

  if(is.null(x$ddf)){
    cat("No detection function\n")
  }else{
    if(as.character(x$formula)[2]!="presence"){
      cat("\nDetection function : ",ddf.model.description(x$ddf),"\n")
    }
  }
  cat("\nFormula: ",as.character(x$formula)[2],"~",
                    as.character(x$formula)[-c(1,2)],"\n")
  cat("\n")

  if(gamm){
    cat("Fitting engine: gamm\n")
  }else{
    cat("Fitting engine: gam\n")
  }

  cat("\n")

  invisible()
}
