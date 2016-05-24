##' Summarizes the gSEM principle 2 result.
##' 
##' summary.sgSEMp1 gives a brief summary about the gSEM Principle 2 analysis. 
##'
##' @title Summary of Semi-gSEM
##' @param object An object of class "sgSEMp2", the returned list from sgSEMp2(). 
##' @param ... Other arguments. Currently not used. 
##' @return NULL. A summary of data and fitting result is printed on screen. 
##'
##' @export
##' 
##' @examples
##' data(acrylic)
##' ans <- sgSEMp2(acrylic)
##' summary(ans)

summary.sgSEMp2 <- function(object, ...){
  cat("gSEM Done \n")
  cat("\n")
  cat("Main predictor:", colnames(object$res.best)[1], "\n")
  cat("Response:", colnames(object$res.best)[2], "\n")
  cat("Intermediate variables:", colnames(object$res.best)[-c(1,2)], "\n")
  cat("\n")
  cat("Chosen models:\n")
  sapply(1:nrow(object$res.print), function(i) {
    paste(object$res.print[i,2], "--->", object$res.print[i,1], "|", "Model:", object$res.print[i,3], 
          "(Adjusted R-square = ", round(as.numeric(object$res.print[i,5]), 5), ")", collapse="\n")
  }
  )
}
