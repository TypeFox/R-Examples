#' \code{initialize} ANFIS object constructor
#'
#' Create the ANFIS object architecture for the trainingSet (X,Y) with full 
#' rules.
#' 
#' @param .Object ANFIS class
#' @param X input matrix with ncol=#inputs and nrow=#individuals 
#' @param Y output matrix with ncol=#output and nrow=#individuals 
#' @param membershipFunction list with the MembershipFunction for each input
#'
#' @return ANFIS object
#'
#' @include Anfis.R
#' @exportMethod initialize
#' @docType methods
#' @name initialize
#' @rdname ANFIS-initialize
#' @aliases initialize,ANFIS-method
#' @seealso \code{\link{ANFIS-class}}
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="initialize",
  signature="ANFIS",
  definition=function(.Object, X, Y, membershipFunction){
  ##Set the different slots
  .Object@X <- X
  .Object@Y <- Y
  .Object@premises <- membershipFunction
  .Object@trainingType <- "Not trained yet"

  ## Generate the complete set of rules, i.e. for 2 input with 2x3 MF=6 rules
  MF<-unlist(lapply(.Object@premises,length))
  .Object@rules<-as.matrix(expand.grid(lapply(MF,function(input){1:input})))
  .Object@rules<-.Object@rules[order(.Object@rules[,1]),,drop=FALSE]

  ## Set the Initial consequents
  .Object@consequents <- matrix(rep(0,ncol(Y)*nrow(.Object@rules)*(ncol(.Object@X)+1)),
        ncol=ncol(Y),nrow=nrow(.Object@rules)*(ncol(.Object@X)+1))
  .Object@errors <- numeric(length=0)

  ##Check the object's validity
  validObject(.Object)

  return(.Object)
})
