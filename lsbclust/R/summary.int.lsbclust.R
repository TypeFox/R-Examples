#' Summary Method for Class "int.lsbclust"
#' 
#' Some goodness-of-fit diagnostics are provided for all three margins.
#' @param object An object of class 'int.lsbclust'.
#' @param digits The number of digits in the printed output.
#' @param \dots Unimplemented.
#' @method summary int.lsbclust
#' @export
summary.int.lsbclust <- function(object, digits = 3, ...){
  
  cat("\tDiagnostics for Interaction Clustering and Biplots\n")
  cat("\n", object$nclust, "clusters;", object$N, "observations\n")
  cat("\nCluster sizes:\n")
  print(table(object$cluster))
  
  ## Overall fit in the SVD
  if(object$fixed == "none") {
    cat("\nVariation accounted for per dimension (rows) for each cluster (columns):\n")
    print(round(object$ofit, digits = min(getOption("digits"), digits)))
  } else {
    cat("\nVariation accounted for per dimension across all clusters:\n")
    print(round(object$ofit, digits = min(getOption("digits"), digits)))
  }
  
  ## Print row fits
  if(object$fixed == "rows") {
    cat("\nVariation accounted for per row across all clusters (sample predictivities):\n")
    print(round(object$rfit, digits = min(getOption("digits"), digits)))
  } else {
    cat("\nVariation accounted for per row per cluster:\n")
    print(round(object$rfit, digits = min(getOption("digits"), digits)))
  }
  
  ## Print column fits
  if(object$fixed == "columns") {
    cat("\nVariation accounted for per column across all clusters (axis predictivities):\n")
    print(round(object$cfit, digits = min(getOption("digits"), digits)))
  }
  else {
    cat("\nVariation accounted for per column per clusters:\n")
    print(round(object$cfit, digits = min(getOption("digits"), digits)))
  }
  
  ## Print loss contributions per person
  cat("\nContributions to total loss per observation (percentage):\n")
  print(summary(object$losscomps*100), digits = min(getOption("digits"), digits))
  
  ## Loss contributions per cluster
  lctb <- with(object, tapply(losscomps, INDEX = cluster, FUN = sum))
  names(lctb) <- seq_len(object$nclust)
  lctb <- c(lctb, Total = object$minloss)
  cat("\nLoss decomposed per cluster:\n")
  print(round(lctb, digits = min(getOption("digits"), digits)))
  
  ## Random starts
  cat("\nObtained loss across", length(object$allloss), "random starts:\n")
  print(summary(object$allloss), digits = min(getOption("digits"), digits))
  
  ## Cluster consistency across random starts
  cat("\nCluster agreement with all other random starts (method = \"", 
      ifelse(is.null(object$call$method) || object$call$method == as.name("method"), formals(int.lsbclust)$method, 
             object$call$method), "\"):\n", sep = "")
  print(summary(drop(object$cl_agreement)), digits = min(getOption("digits"), digits))
  cat("\nSee ?clue::cl_agreement for the interpretation of 'method'.")
}