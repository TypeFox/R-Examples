###################################################################################
##                                MixmodDAResults.R                              ##
###################################################################################

###################################################################################
##' @include global.R
##' @include MixmodResults.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{MixmodDAResults}}] class
##'
##' This is a class to contain results after a discriminant analysis with MIXMOD. Inherits the [\code{\linkS4class{MixmodResults}}] class.
##'  
##' \describe{
##'   \item{CVLabel}{vector of integers containing labels defined by cross validation.}
##'   \item{CVClassification}{classification table after cross validation.}
##'   \item{MAPErrorRate}{error rate done by MAP algorithm.}
##'   \item{MAPClassification}{classification table after MAP algorithm.}
##' }
##'
##' @examples
##'   getSlots("MixmodDAResults")
##'
##' @name MixmodDAResults-class
##' @rdname MixmodDAResults-class
##' @exportClass MixmodDAResults
##'
setClass(
  Class="MixmodDAResults",
  representation=representation(
    CVLabel = "integer",
    CVClassification = "matrix",
    MAPErrorRate = "numeric",
    MAPClassification = "matrix"
  ),
  contains = c("MixmodResults"),
  prototype=prototype(
    CVLabel = integer(0),
    MAPErrorRate = numeric(0),
    CVClassification = matrix(nrow=0,ncol=0),
    MAPClassification = matrix(nrow=0,ncol=0)
  )
)
###################################################################################

###################################################################################
##' @rdname print-methods
##' @aliases print print,MixmodDAResults-method
##'
setMethod(
  f="print",
  signature=c("MixmodDAResults"),
  function(x,...){ 
    callNextMethod()
    if ( length(x@CVLabel)>0){
      cat("* Classification with CV:\n")
      cat("           |");
      line=c("-----------")
      for( i in 1:x@nbCluster ){
        cat(" Cluster", i, "|")
        line=c(line,"-----------")
      }
      cat("\n");cat(line,"\n")
      for( i in 1:x@nbCluster ){
        cat(" Cluster", i, "|")
        for( j in 1:x@nbCluster ){
          cat(formatC(x@CVClassification[i,j],digits=9,format="d"),"|")
        }
        cat("\n")
      }
      cat(line,"\n")
      cat("* Error rate with CV = ", formatC((1-x@criterionValue)*100,digits=2,format="f"),"%\n\n")
    }
    cat("* Classification with MAP:\n")
    cat("           |");
    line=c("-----------")
    for( i in 1:x@nbCluster ){
      cat(" Cluster", i, "|")
      line=c(line,"-----------")
    }
    cat("\n");cat(line,"\n")
    for( i in 1:x@nbCluster ){
      cat(" Cluster", i, "|")
      for( j in 1:x@nbCluster ){
        cat(formatC(x@MAPClassification[i,j],digits=9,format="d"),"|")
      }
      cat("\n")
    }
    cat(line,"\n")
    cat("* Error rate with MAP = ", formatC(x@MAPErrorRate*100,digits=2,format="f"),"%\n")

    cat("****************************************\n")
  }
)
###################################################################################

###################################################################################
##' @rdname show-methods
##' @aliases show show,MixmodDAResults-method
##'
setMethod(
  f="show",
  signature=c("MixmodDAResults"),
  function(object){
    callNextMethod()
    if ( length(object@CVLabel)>0){
      cat("* Classification with CV:\n")
      cat("           |");
      line=c("-----------")
      for( i in 1:object@nbCluster ){
        cat(" Cluster", i, "|")
        line=c(line,"-----------")
      }
      cat("\n");cat(line,"\n")
      for( i in 1:object@nbCluster ){
        cat(" Cluster", i, "|")
        for( j in 1:object@nbCluster ){
          cat(formatC(object@CVClassification[i,j],digits=9,format="d"),"|")
        }
        cat("\n")
      }
      cat(line,"\n")
      cat("* Error rate with CV = ", formatC((1-object@criterionValue)*100,digits=2,format="f"),"%\n\n")
    }
    cat("* Classification with MAP:\n")
    cat("           |");
    line=c("-----------")
    for( i in 1:object@nbCluster ){
      cat(" Cluster", i, "|")
      line=c(line,"-----------")
    }
    cat("\n");cat(line,"\n")
    for( i in 1:object@nbCluster ){
      cat(" Cluster", i, "|")
      for( j in 1:object@nbCluster ){
        cat(formatC(object@MAPClassification[i,j],digits=9,format="d"),"|")
      }
      cat("\n")
    }
    cat(line,"\n")
    cat("* Error rate with MAP = ", formatC(object@MAPErrorRate*100,digits=2,format="f"),"%\n")

    cat("****************************************\n")
  }
)
###################################################################################




