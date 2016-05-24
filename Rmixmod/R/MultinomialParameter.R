###################################################################################
##                        MultinomialParameter.R                                 ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Parameter.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{MultinomialParameter}}] class
##' 
##' This class defines parameters of a Multinomial Mixture Model. Inherits the [\code{\linkS4class{Parameter}}] class.
##' 
##' \describe{
##'   \item{center}{a numeric vector containing center of each cluster.}
##'   \item{scatter}{a vector of matrix containing dispersion matrix of each cluster.}
##'   \item{factor}{a character vector containing the modalities.}
##' }
##'
##' @examples
##'   new("MultinomialParameter")
##'
##'   getSlots("MultinomialParameter")
##' 
##' @name MultinomialParameter-class
##' @rdname MultinomialParameter-class
##' @exportClass MultinomialParameter
##'
setClass(
    Class="MultinomialParameter",
    representation=representation(
        center = "numeric",
        scatter = "numeric",
        factor = "factor"
    ),
    prototype=prototype(
        center = numeric(0),
        scatter =  numeric(0),
        factor = factor()
    ),
    contains=c("Parameter")
)
###################################################################################



###################################################################################
##' @rdname print-methods
##' @aliases print print,MultinomialParameter-method
##'
setMethod(
    f="print",
    signature=c("MultinomialParameter"),
    function(x,...){
      if(length(x@proportions)>0){
        cat("****************************************\n")
        cat("* number of modalities   = ", x@factor, "\n")
        for(i in 1:length(x@proportions)){
          cat("*** Cluster",i,"\n")
          cat("* proportion = ", formatC(x@proportions[i],digits=4,format="f"), "\n")
          cat("* center     = ", formatC(x@center[i,],digits=4,format="f"), "\n")
          cat("* scatter    = |",formatC(x@scatter[[i]][1,1:x@factor[1]],digits=4,width=10,format="f"),"|\n")
          for ( j in 2:nrow(x@scatter[[i]])){
            cat("               |", formatC(x@scatter[[i]][j,1:x@factor[j]],digits=4,width=10,format="f"),"|\n")
          }
        }
        cat("****************************************\n")
      }
    }
)
###################################################################################


###################################################################################
##' @rdname show-methods
##' @aliases show show,MultinomialParameter-method
##'
setMethod(
    f="show",
    signature=c("MultinomialParameter"),
    function(object){
      if(length(object@proportions)>0){
        cat("****************************************\n")
        cat("* number of modalities = ", object@factor, "\n")
        for(i in 1:length(object@proportions)){
          cat("*** Cluster",i,"\n")
          cat("* proportion = ", formatC(object@proportions[i],digits=4,format="f"), "\n")
          cat("* center     = ", formatC(object@center[i,],digits=4,format="f"), "\n")
          if ( length(object@scatter) > 1 ){
            if( length(dim(object@scatter)) == 2 ){
              cat("* scatter    = ", formatC(object@scatter[i,],digits=4,format="f"), "\n")
            }else if( length(object@scatter[[i]]) > 1 ){
              if ( nrow(object@scatter[[i]])>1 ){
                cat("* scatter    = |",formatC(object@scatter[[i]][1,1:object@factor[1]],digits=4,width=10,format="f"),"|\n")
                for ( j in 2:nrow(object@scatter[[i]])){
                  cat("               |", formatC(object@scatter[[i]][j,1:object@factor[j]],digits=4,width=10,format="f"),"|\n")
                }
              }else{
                cat("* scatter    = ",formatC(object@scatter[[i]],digits=4,format="f"),"\n")
              }
            }else{
              cat("* scatter    = ", formatC(object@scatter[i],digits=4,format="f"), "\n")
            }
          }else{
            cat("* scatter    = ", formatC(object@scatter,digits=4,format="f")," \n")
          }
        }
        cat("****************************************\n")
      }
    }
)
###################################################################################


###################################################################################
##' @rdname summary-methods
##' @aliases summary summary,MultinomialParameter-method
##'
setMethod(
  f="summary",
  signature=c("MultinomialParameter"),
  function(object, ...){
    if(length(object@proportions)>0){
      for ( k in 1:length(object@proportions) ){
        cat("*                  Cluster ", k,": \n")
        cat("                         Proportion = ",formatC(object@proportions[k],digits=4,format="f"),"\n")
        cat("                             Center = ",formatC(object@center[k,],digits=4,format="f"),"\n")
        if ( length(object@scatter) > 1 ){
          if( length(dim(object@scatter)) == 2 ){
            cat("                            Scatter = ", formatC(object@scatter[k,],digits=4,format="f"), "\n")
          }else if( length(object@scatter[[k]]) > 1 ){
            if ( nrow(object@scatter[[k]])>1 ){
              cat("                            Scatter = |",formatC(object@scatter[[k]][1,1:object@factor[1]],digits=4,width=10,format="f"),"|\n")
              for ( j in 2:nrow(object@scatter[[k]])){
                cat("                                      |", formatC(object@scatter[[k]][j,1:object@factor[j]],digits=4,width=10,format="f"),"|\n")
              }
            }else{
              cat("                            Scatter = ",formatC(object@scatter[[k]],digits=4,format="f"),"\n")
            }
          }else{
            cat("                            Scatter = ", formatC(object@scatter[k],digits=4,format="f"), "\n")
          }
        }else{
          cat("                            Scatter = ", formatC(object@scatter,digits=4,format="f")," \n")
        }
      }
    }
  }
)
###################################################################################


###################################################################################
##' @rdname extract-methods
##' @aliases [,MultinomialParameter-method
##'
setMethod(
  f="[", 
  signature(x = "MultinomialParameter"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "proportions"={return(x@proportions)},
        "center"={return(x@center)},
        "scatter"={return(x@scatter)},
        stop("This attribute doesn't exist !")
      )
    }else{
        switch(EXPR=i,
        "proportions"={return(x@proportions[j])},
        "center"={return(x@center[j,])},
        "scatter"={return(x@scatter[j])},
        stop("This attribute doesn't exist !")
        )
    }
  }
)
###################################################################################


###################################################################################
##' @name [
##' @rdname extract-methods
##' @aliases [<-,MultinomialParameter-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "MultinomialParameter"), 
  definition=function(x,i,j,value){
    if ( missing(j) ){
      switch(EXPR=i,
      "proportions"={x@proportions<-value},
      "center"={x@center<-value},
      "scatter"={x@scatter<-value},
      stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
      "proportions"={x@proportions[j]<-value},
      "center"={x@center[j,]<-value},
      "scatter"={x@scatter[j]<-value},
      stop("This attribute doesn't exist !")
      )
    }
    validObject(x)
    return(x)
  }
)
###################################################################################

