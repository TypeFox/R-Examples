#' \code{Print} and \code{Show} an ANFIS object
#'
#' Generic Print/Show Method for ANFIS class output visualization. 
#'
#' @param x ANFIS class object
#' @param object ANFIS class object
#' @param ... not used but included for generic print compatibility
#'
#' @return console output of the object
#'
#' @include Anfis-getters.R
#' @exportMethod print
#' @docType methods
#' @name print
#' @rdname ANFIS-printshow
#' @aliases print,ANFIS-method
#' @seealso \code{\link{ANFIS-class}}
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="print", signature="ANFIS", definition=function(x, ...){
  show(x)
  cat("X: \n"); print(x@X[1:2,,drop=FALSE])
  cat("Y: \n"); print(x@Y[1:2,,drop=FALSE])
  cat("Premises, first 2 per input: \n") 
  print(lapply(x@premises,function(input){ lapply(input[1:2],"[") }))
  cat("Rules, first 2: \n"); print(x@rules[1:2,,drop=FALSE])
  cat("Consequents, first 2: \n"); print(x@consequents[1:2,,drop=FALSE])
})
#' @exportMethod show
#' @docType methods
#' @name show
#' @rdname ANFIS-printshow
#' @inheritParams print
#' @aliases show,ANFIS-method
setMethod(f="show", signature="ANFIS", definition=function(object){
  ##Print class header
  cat("ANFIS network \n")
  cat("Trainning Set: \n")
  cat("\t dim(x)=", gsub(", ","x",toString(dim(object@X))),"\n") 
  cat("\t dim(y)=",gsub(", ","x",toString(dim(object@Y))),"\n")
  cat("Arquitecture: ", ncol(object@X), "(",
    gsub(", ","x", toString(unlist(lapply(object@premises,length)))),
    ") -", nrow(object@rules), "-", length(object@consequents), "(", 
    gsub(", ", "x", toString(dim(object@consequents))), ") -", 
    ncol(object@consequents),"\n")
  
  ##Check if network is trained
  if(length(object@errors)==0){
    cat("Network not trained yet \n")
  }else{
    cat("Last training error: ", object@errors[length(object@errors)],"\n")
  }
})
