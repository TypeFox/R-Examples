###################################################################################
##                             CompositeParameter.R                               ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Parameter.R
##' @include GaussianParameter.R
##' @include MultinomialParameter.R

NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{CompositeParameter}}] class
##' 
##' This class defines parameters of a Heterogeneous Mixture Model. Inherits the [\code{\linkS4class{Parameter}}] class.
##' 
##' \describe{
##'   \item{g_parameter}{an object of class [\code{\linkS4class{GaussianParameter}}] }
##'   \item{m_parameter}{an object of class [\code{\linkS4class{MultinomialParameter}}] }
##'   \item{factor}{numeric for factor}
##' }
##'
##' @examples
##'   new("CompositeParameter")
##'
##'   getSlots("CompositeParameter")
##' 
##' @name CompositeParameter-class
##' @rdname CompositeParameter-class
##' @exportClass CompositeParameter
##'
setClass(
  Class="CompositeParameter",
  representation=representation(
    g_parameter = "GaussianParameter",
    m_parameter = "MultinomialParameter",
    factor = "numeric"
  ),
  contains=c("Parameter")
)
###################################################################################


###################################################################################
##' @rdname print-methods
##' @aliases print print,CompositeParameter-method
##'
setMethod(
  f="print",
  signature=c("CompositeParameter"),
  function(x,...){
    if(length(x@proportions)>0){
      cat("Gaussian Parameters\n")
      print(x@g_parameter)
      cat("Multinomial Parameters\n")
      print(x@m_parameter)
    }
  }
)
###################################################################################


###################################################################################
##' @rdname show-methods
##' @aliases show show,CompositeParameter-method
##'
setMethod(
  f="show",
  signature=c("CompositeParameter"),
  function(object){
    if(length(object@proportions)>0){
      cat("Gaussian Parameters\n")
      show(object@g_parameter)
      cat("Multinomial Parameters\n")
      show(object@m_parameter)
    }
  }
)
###################################################################################


###################################################################################
##' @rdname summary-methods
##' @aliases summary summary,CompositeParameter-method
##'
setMethod(
  f="summary",
  signature=c("CompositeParameter"),
  function(object, ...){
    if(length(object@proportions)>0){
      cat("Gaussian Parameters\n")
      summary(object@g_parameter)
      cat("Multinomial Parameters\n")
      summary(object@m_parameter)
    }
  }
)
###################################################################################


###################################################################################
##' @rdname extract-methods
##' @aliases [,CompositeParameter-method
##'
setMethod(
  f="[", 
  signature(x = "CompositeParameter"),
  definition=function(x,i,j,drop){
      switch(EXPR=i,
             "g_parameter"={return(x@g_parameter)},
             "m_parameter"={return(x@m_parameter)},
             stop("This attribute doesn't exist !")
             )
  }
)
###################################################################################



###################################################################################
##' @name [
##' @rdname extract-methods
##' @aliases [<-,CompositeParameter-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "CompositeParameter"), 
  definition=function(x,i,j,value){
      switch(EXPR=i,
             "g_parameter"={return(x@g_parameter)},
             "m_parameter"={return(x@m_parameter)},
             stop("This attribute doesn't exist !")
      )
    validObject(x)
    return(x)
  }
)
###################################################################################
