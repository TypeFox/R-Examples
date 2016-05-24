#' Automatic Dummy Variable Creation with Support for Predictive Contexts
#'
#' \code{dummy} creates dummy variables of all the factors and character vectors in a data frame. It also supports settings in which the user only wants to compute dummies for the categorical values that were present in another data set. This is especially useful in the context of predictive modeling, in which the new (test) data has more or other categories than the training data.
#'
#' @param x a data frame containing at least one factor or character vector
#' @param p Only relevant if object is NULL. Select the top p values in terms of frequency. Either "all" (all categories in all variables), an integer scalar (top p categories in all variables), or a vector of integers (number of top categories per variable in order of appearance).
#' @param object output of the \code{categories} function. This parameter is to be used when dummies should be created only of categories present in another data set (e.g., training set)
#' @param int should the dummies be integers (TRUE) or factors (FALSE)
#' @param verbose logical. Used to show progress
#' @examples
#' #create toy data
#' (traindata <- data.frame(var1=as.factor(c("a","b","b","c")),
#'                          var2=as.factor(c(1,1,2,3)),
#'                          var3=c("val1","val2","val3","val3"),
#'                          stringsAsFactors=FALSE))
#' (newdata <- data.frame(var1=as.factor(c("a","b","b","c","d","d")),
#'                        var2=as.factor(c(1,1,2,3,4,5)),
#'                        var3=c("val1","val2","val3","val3","val4","val4"),
#'                        stringsAsFactors=FALSE))
#' #create dummies of training set
#' (dummies_train <- dummy(x=traindata))
#' #create dummies of new set
#' (dummies_new <- dummy(x=newdata))
#'
#' #how many new dummy variables should not have been created?
#' sum(! colnames(dummies_new) %in% colnames(dummies_train))
#'
#' #create dummies of new set using categories found in training set
#' (dummies_new <- dummy(x=newdata,object=categories(traindata,p="all")))
#'
#' #how many new dummy variables should not have be created?
#' sum(! colnames(dummies_new) %in% colnames(dummies_train))
#'
#'
#' #create dummies of training set, 
#' #using the top 2 categories of all variables found in the training data
#' dummy(x=traindata,p=2)
#' 
#' #create dummies of training set, 
#' #using respectively the top 2,3 and 1 categories of the three 
#' #variables found in training data
#' dummy(x=traindata,p=c(2,3,1))
#' 
#' #create all dummies of training data
#' dummy(x=traindata)
#' 
#'
#' @seealso \code{\link{categories}}
#' @return  A data frame containing dummy variables
#' @author Authors: Michel Ballings, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
dummy <- function(x,p="all", object=NULL,int=FALSE,verbose=FALSE){
  colnames(x) <- make.names(colnames(x),TRUE)
  if(is.null(object)) object <- categories(x,p=p)
  ans <- list()
  #for each variable
  len <- length(object)
  for (i in 1:len){
    if (verbose) cat(round((i*100)/len,0),"% of variables processed \n")
    #for each value an ifelse
    ans[[i]] <- data.frame(sapply(object[[i]], function(z) {
      if (int==FALSE) {
        z <- as.factor(ifelse(x[,names(object)[i]] == z,1,0))
      } else if (int==TRUE) {
        z <- as.integer(ifelse(x[,names(object)[i]] == z,1,0))
      }
      z}
    )                    
    )
    colnames(ans[[i]]) <- make.names(paste(names(object)[i],object[[i]],sep="_"),TRUE)
    
  }
  do.call(cbind,ans)
}