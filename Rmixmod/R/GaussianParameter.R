###################################################################################
##                             GaussianParameter.R                               ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Parameter.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{GaussianParameter}}] class
##' 
##' This class defines parameters of a Gaussian Mixture Model. Inherits the [\code{\linkS4class{Parameter}}] class.
##' 
##' \describe{
##'   \item{mean}{a matrix containing mean of each cluster.}
##'   \item{variance}{a list of matrices containing variance matrix of each cluster.}
##' }
##'
##' @examples
##'   new("GaussianParameter")
##'
##'   getSlots("GaussianParameter")
##' 
##' @name GaussianParameter-class
##' @rdname GaussianParameter-class
##' @exportClass GaussianParameter
##'
setClass(
    Class="GaussianParameter",
    representation=representation(
        mean = "matrix",
        variance = "list"
    ),
    prototype=prototype(
        mean = matrix(0),
        variance = list(0)
    ),
    contains=c("Parameter")
)
###################################################################################


###################################################################################
##' @rdname print-methods
##' @aliases print print,GaussianParameter-method
##'
setMethod(
    f="print",
    signature=c("GaussianParameter"),
    function(x,...){
      if(length(x@proportions)>0){
        cat("****************************************\n")
        for(k in 1:length(x@proportions)){
          cat("*** Cluster",k,"\n")
          cat("* proportion = ", formatC(x@proportions[k],digits=4,format="f"), "\n")
          cat("* means      = ", formatC(x@mean[k,],digits=4,format="f"), "\n")
          if ( nrow(x@variance[[k]])>1 ){
            cat("* variances  = |",formatC(x@variance[[k]][1,],digits=4,width=10,format="f"),"|\n")
            for ( i in 2:nrow(x@variance[[k]])){
              cat("               |", formatC(x@variance[[k]][i,],digits=4,width=10,format="f"),"|\n")
            }
          }else{
            cat("* variances  = ",formatC(x@variance[[k]],digits=4,format="f"),"\n")
          }
        }
        cat("****************************************\n")
      }
    }
)
###################################################################################


###################################################################################
##' @rdname show-methods
##' @aliases show show,GaussianParameter-method
##'
setMethod(
    f="show",
    signature=c("GaussianParameter"),
    function(object){
      if(length(object@proportions)>0){
        cat("****************************************\n")
        for(k in 1:length(object@proportions)){
          cat("*** Cluster",k,"\n")
          cat("* proportion = ", formatC(object@proportions[k],digits=4,format="f"), "\n")
          cat("* means      = ", formatC(object@mean[k,],digits=4,format="f"), "\n")
          if ( nrow(object@variance[[k]])>1 ){
            cat("* variances  = |",formatC(object@variance[[k]][1,],digits=4,width=10,format="f"),"|\n")
            for ( i in 2:nrow(object@variance[[k]])){
              cat("               |", formatC(object@variance[[k]][i,],digits=4,width=10,format="f"),"|\n")
            }
          }else{
            cat("* variances  = ",formatC(object@variance[[k]],digits=4,format="f"),"\n")
          }
        }
        cat("****************************************\n")
      }
    }
)
###################################################################################


###################################################################################
##' @rdname summary-methods
##' @aliases summary summary,GaussianParameter-method
##'
setMethod(
  f="summary",
  signature=c("GaussianParameter"),
  function(object, ...){
    if(length(object@proportions)>0){
      for ( k in 1:length(object@proportions) ){
        cat("*                  Cluster ", k,": \n")
        cat("                         Proportion = ",formatC(object@proportions[k],digits=4,format="f"),"\n")
        cat("                              Means = ",formatC(object@mean[k,],digits=4,format="f"),"\n")
        if ( nrow(object@variance[[k]])>1 ){
          cat("                          Variances = |",formatC(object@variance[[k]][1,],digits=4,width=10,format="f"),"|\n")
          for ( i in 2:nrow(object@variance[[k]])){
            cat("                                      |", formatC(object@variance[[k]][i,],digits=4,width=10,format="f"),"|\n")
          }
        }else{
          cat("                          Variances = ",formatC(object@variance[[k]],digits=4,format="f"),"\n")
        }
      }
    }
  }
)
###################################################################################


###################################################################################
##' @rdname extract-methods
##' @aliases [,GaussianParameter-method
##'
setMethod(
  f="[", 
  signature(x = "GaussianParameter"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "proportions"={return(x@proportions)},
        "mean"={return(x@mean)},
        "variance"={return(x@variance)},
        stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
        "proportions"={return(x@proportions[j])},
        "mean"={return(x@mean[j,])},
        "variance"={return(x@variance[[j]])},
        stop("This attribute doesn't exist !")
      )
    }
  }
)
###################################################################################



###################################################################################
##' @name [
##' @rdname extract-methods
##' @aliases [<-,GaussianParameter-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "GaussianParameter"), 
  definition=function(x,i,j,value){
    if ( missing(j) ){
      switch(EXPR=i,
        "proportions"={x@proportions<-value},
        "mean"={x@mean<-value},
        "variance"={x@variance<-value},
        stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
        "proportions"={x@proportions[j]<-value},
        "mean"={x@mean[j,]<-value},
        "variance"={x@variance[[j]]<-value},
        stop("This attribute doesn't exist !")
        )
    }
    validObject(x)
    return(x)
  }
)
###################################################################################

