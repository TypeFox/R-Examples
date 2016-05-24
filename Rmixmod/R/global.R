###################################################################################
##                                global.R                                       ##
###################################################################################

###################################################################################
##' Define function to check an integer
##'
##' @param x a numeric
##' @param tol a real which defines the tolerance
##'
##' @return a logical. TRUE if \code{x} is an integer, FALSE otherwise.
##' 
##' @keywords internal
##'
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
###################################################################################


###################################################################################
##' Define function to transform a vector of modalities into matrix of binary data
##'
##' @param x a numeric
##'
##' @return a matrix of binary data.
##' 
##' @keywords internal
##'
vector2binary <- function (x) {
  x<-as.factor(x)
  y <- matrix(0, nrow=length(x), ncol=length(levels(x)))
  for (i in 1:length(x)) {
    y[i, as.numeric(x[i])] <- 1
  }
  return(y)
}
###################################################################################

###################################################################################
##' Define function to transform a matrix of modalities into matrix of binary data
##'
##' @param x a matrix
##'
##' @return a matrix of binary data.
##' 
##' @keywords internal
##'
matrix2binary <- function (x) {
  if( is.vector(x) )
  y <- vector2binary(x)
  else {
    y <- NULL
    y.names <- NULL
    for (i in 1:length(x)) {
      y <- cbind(y, vector2binary(x[,i]))
      y.names <- c( y.names, paste(names(x)[i], levels(as.factor(x[,i])), sep='') )
    }
  }
  colnames(y) <- y.names
  return(y)
}
###################################################################################

###################################################################################
##' Say if a data frame is quantitative, qualitative or composite
##'
##' @param x a vector, a factor or a data frame
##'
##' @return a string with the data type
##' 
##' @export
##'
is.dataType <- function (x) {
  # loof after missing values
  if ( sum(is.na(x)) ) stop("data set contains missing values")
  if ( is.vector(x) | is.factor(x) ){
    if ( is.double(x) ) return("quantitative")
    if ( nlevels(as.factor(x))/length(x)>.25 ) 
      warning("more than 25% of the observations have different modalities")
    return("qualitative")
  }
  else{
    if (sum(sapply(x,is.numeric))==ncol(x))
      return("quantitative")
    if (sum(sapply(x,is.factor))==ncol(x))
      return("qualitative")
    else
      return ("composite")
  }
}
###################################################################################

###################################################################################
##' Say if a data frame contains only qualitative variables.
##'
##' @param x a vector or a matrix or a data frame
##'
##' @return a boolean
##' 
##' @export
##'
isQualitative <- function (x) {
  # loof after missing values
  if ( sum(is.na(x)) ) stop("data set contains missing values")
  if ( is.vector(x) | is.factor(x) ){
    if ( is.double(x) ) return(FALSE)
    if ( nlevels(as.factor(x))/length(x)>.25 ) 
      warning("more than 25% of the observations have different modalities")
    return(TRUE)
  }
  else{
    # loop over columns to check whether type is factor
    for ( j in 1:ncol(x) ){
     if ( is.double(x[,j]) ) {return(FALSE)}
     if ( nlevels(as.factor(x[,j]))/length(x[,j])>.25 ) 
       warning(paste("more than 25% of the observations have different modalities for variable named",names(x)[j]))
    }
  }
  return(TRUE)
}
###################################################################################

###################################################################################
##' Convert a data frame containing integers to a qualitative data set with factors.
##'
##' @param x a vector or a matrix or a data frame
##'
##' @return a qualitative data set with factors
##' 
##' @export
##'
asQualitative <- function (x) {
  if ( isQualitative(x) ){
    if ( is.vector(x) ) x<-as.factor(x)
    else if ( is.data.frame(x) | is.matrix(x)  ){
      # loop over columns to check whether type is factor
      for ( j in 1: ncol(x) ){ x[,j] <- as.factor(x[,j]) }
    }
  }
  return(x)
}
###################################################################################


###################################################################################
##' Get the number of modalities for each column of a categorical data set
##'
##' @param x a vector or a matrix or a data frame containing factors
##'
##' @return a vector containing the number of modalities for each column
##' 
##' @export
##'
nbFactorFromData <- function (x) {
  # set x as a qualitative data set
  x<-asQualitative(x)
  if ( is.vector(x) | is.factor(x) ){
    return(nlevels(x))
  }else{
    return(sapply(x,nlevels))
  }
}

###################################################################################


###################################################################################
##' Sorting results of a [\code{\linkS4class{Mixmod}}] object by a given criterion
##'
##' After calling the mixmodCluster() or mixmodLearn() method, results will be sorted
##' into ascending order according to the first given criterion (descending order for CV criterion). 
##' This method is able to reorder the list of results according to a given criterion.
##'
##' @param object a [\code{\linkS4class{Mixmod}}] object
##' @param criterion a string containing the criterion name
##'
##' @return a modified [\code{\linkS4class{Mixmod}}] object
##' 
##' @exportMethod sortByCriterion
##' @docType methods
##' @name sortByCriterion
##' @rdname sortByCriterion-methods
##' @examples
##'   x<-mixmodCluster(iris[1:4],2:10,criterion=c("BIC","ICL"))
##'   icl<-sortByCriterion(x,"ICL")
##'   icl["results"]
##'
setGeneric(
  name="sortByCriterion",
  def=function(object,criterion){
    standardGeneric("sortByCriterion")
  }
)
###################################################################################


###################################################################################
##' Print a Rmixmod class to standard output.
##'
##' @param x a Rmixmod object: a \code{\linkS4class{Strategy}}, a \code{\linkS4class{Model}}, a \code{\linkS4class{GaussianParameter}}, a \code{\linkS4class{MultinomialParameter}}, a \code{\linkS4class{MixmodResults}}, a \code{\linkS4class{MixmodCluster}}, a \code{\linkS4class{MixmodLearn}} or a \code{\linkS4class{MixmodPredict}}.
##' @param ... further arguments passed to or from other methods
##'
##' @return NULL. Prints to standard out.
##'
##' @name print
##' @rdname print-methods
##' @docType methods
##' @exportMethod print
##'
##' @seealso \code{\link{print}}
##' @examples
##'   ## for strategy
##'   strategy <- mixmodStrategy()
##'   print(strategy)
##'
##'   ## for Gaussian models
##'   gmodel <- mixmodGaussianModel()
##'   print(gmodel)
##'   ## for multinomial models
##'   mmodel <- mixmodMultinomialModel()
##'   print(mmodel)
##'
##'   ## for clustering
##'   data(geyser)
##'   xem <- mixmodCluster(geyser,3)
##'   print(xem)
##'   ## for Gaussian parameters
##'   print(xem["bestResult"]["parameters"])
##'
##'   ## for discriminant analysis
##'   # start by extract 10 observations from iris data set
##'   iris.partition<-sample(1:nrow(iris),10)
##'   # then run a mixmodLearn() analysis without those 10 observations
##'   learn<-mixmodLearn(iris[-iris.partition,1:4], iris$Species[-iris.partition])
##'   # print learn results
##'   print(learn)
##'   # create a MixmodPredict to predict those 10 observations
##'   prediction <- mixmodPredict(data=iris[iris.partition,1:4], classificationRule=learn["bestResult"])
##'   # print prediction results
##'   print(prediction)
##'
NULL

###################################################################################
##' Show description of a Rmixmod class to standard output.
##'
##' @param object a Rmixmod object: a \code{\linkS4class{Strategy}}, a \code{\linkS4class{Model}}, a \code{\linkS4class{GaussianParameter}}, a \code{\linkS4class{MultinomialParameter}}, a \code{\linkS4class{MixmodResults}}, a \code{\linkS4class{MixmodCluster}}, a \code{\linkS4class{MixmodLearn}} or a \code{\linkS4class{MixmodPredict}}.
##'
##' @return NULL. Prints to standard out.
##'
##' @importFrom methods show
##' @name show
##' @docType methods
##' @rdname show-methods
##' @exportMethod show
##'
##' @seealso \code{\link{show}}
##' @examples
##'   ## for strategy
##'   strategy <- mixmodStrategy()
##'   show(strategy)
##'
##'   ## for Gaussian models
##'   gmodel <- mixmodGaussianModel()
##'   show(gmodel)
##'   ## for multinomial models
##'   mmodel <- mixmodMultinomialModel()
##'   show(mmodel)
##'
##'   ## for clustering
##'   data(geyser)
##'   xem <- mixmodCluster(geyser,3)
##'   show(xem)
##'   ## for Gaussian parameters
##'   show(xem["bestResult"]["parameters"])
##'
##'   ## for discriminant analysis
##'   # start by extract 10 observations from iris data set
##'   iris.partition<-sample(1:nrow(iris),10)
##'   # then run a mixmodLearn() analysis without those 10 observations
##'   learn<-mixmodLearn(iris[-iris.partition,1:4], iris$Species[-iris.partition])
##'   # create a MixmodPredict to predict those 10 observations
##'   prediction <- mixmodPredict(data=iris[iris.partition,1:4], classificationRule=learn["bestResult"])
##'   # show results
##'   show(prediction)
##'
NULL


###################################################################################
##' Produce result summaries of a Rmixmod class 
##'
##' @param x Either a \code{\linkS4class{Mixmod}} object, a \code{\linkS4class{Parameter}} object or a \code{\linkS4class{MixmodResults}} object.
##' @param ... further arguments passed to or from other methods
##'
##' @return NULL. Summaries to standard out.
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##' @exportMethod summary
##'
##' @seealso \code{\link{summary}}
##' @examples
##'   data(geyser)
##'   xem <- mixmodCluster(geyser,3)
##'   summary(xem)
##'   summary(xem["bestResult"])
##'   summary(xem["bestResult"]["parameters"])
##'
NULL


###################################################################################
##' Extract parts of a Rmixmod class
##' 
##' @param x object from which to extract element(s) or in which to replace element(s).
##' @param i the name of the element we want to extract or replace.
##' @param j if the element designing by i is complex, j specifying elements to extract or replace.
##' @param drop For matrices and arrays.  If TRUE the result is coerced to the lowest possible dimension (see the examples).  This only works for extracting elements, not for the replacement.  See drop for further details.
##'
##' @name [
##' @docType methods
##' @rdname extract-methods
##'
NULL

