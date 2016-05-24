##' Summarizes the gSEM principle 1 result.
##' 
##' summary.sgSEMp1 gives a summary about the gSEm-Principle 1 analysis. 
##'
##' @title Summary of Semi-gSEM
##' @param object An object of class "sgSEMp1", the returned list from sgSEMp1(). 
##' @param ... Other arguments. Currently not used. 
##' @return NULL. A summary of data and fitting result is printed on screen. 
##'
##' @export
##' 
##' @examples
##' data(acrylic)
##' ans <- sgSEMp1(acrylic)
##' summary(ans)

summary.sgSEMp1 <- function(object, ...){
  cat("gSEM principle 1 Done \n")
  cat("\n")
  cat("Main predictor:", colnames(object$bestModels)[1], "\n")
  cat("Response:", colnames(object$bestModels)[2], "\n")
  cat("Intermediate variables:", colnames(object$bestModels)[-c(1,2)], "\n")
  cat("\n")
  cat("Chosen models:\n")
  sapply(1:nrow(object$table), function(i) {
      paste(object$table[i,2], "--->", object$table[i,1],
            "|", object$table[i,3], "(Adjusted R-square = ",
            round(as.numeric(object$table[i,5]), 5), ")", collapse="\n")
  }
  )
}
