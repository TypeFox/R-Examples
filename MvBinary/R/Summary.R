########################################################################################################################
## Surcharge des fonctions summary et print pour les objects de classe S4 MvBinaryResult
########################################################################################################################

#'
#' Summary function.
#' 
#' This function gives the summary of output from \code{MvBinaryEstim}.
#' 
#' @param object output object from \code{\link{MvBinaryEstim}}.
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' 
#' 
NULL

#' @rdname summary-methods
#' @aliases summary summary,MvBinaryResult-method
#'
setMethod(
  f="summary",
  signature = c("MvBinaryResult"),
  definition = function(object){
    cat("****************************************************************************************\n")
    cat("The model contains", max(object@blocks), "blocks\n")
    cat("Its log-likelihood is", object@loglike, "\n")
    cat("Its BIC criterion value is", object@bic, "\n")
    cat("The model requires", object@nbparam, "parameters\n")
    cat("\n")
    cat("****************************************************************************************\n")
    cat("The blocks are defined as follows\n")
    for (k in 1:max(object@blocks))
      cat("Block", k, "contains the following", sum(object@blocks==k), "variables:\n", names(object@blocks)[which(object@blocks==k)],"\n\n")    
    cat("****************************************************************************************")
  }
)



#'
#' Summary function.
#' 
#' This function prints the parameters resulting from \code{MvBinaryEstim}.
#' 
#' @param x output object from \code{\link{MvBinaryEstim}}.
#' 
#' @name print
#' @rdname print-methods
#' @docType methods
#' @exportMethod print
#' 
#' 
NULL

#' @rdname print-methods
#' @aliases print print,MvBinaryResult-method
#'
setMethod(
  f="print",
  signature = c("MvBinaryResult"),
  definition = function(x){
    summary(x)    
    cat("\n")
    for (k in 1:max(x@blocks)){
      tmp <- cbind(x@alpha[which(x@blocks==k)], x@epsilon[which(x@blocks==k)], x@delta[which(x@blocks==k)])
      colnames(tmp) <- c("alpha", "epsilon", "delta")
      rownames(tmp) <- names(x@blocks)[which(x@blocks==k)]
      cat("Parameters of Block",k,"\n")
      print(round(tmp,3))
      cat("\n")
    }
  }
)