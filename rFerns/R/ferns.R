#    R part of rFerns
#
#    Copyright 2011-2015 Miron B. Kursa
#
#    This file is part of rFerns R package.
#
#rFerns is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#rFerns is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with rFerns. If not, see http://www.gnu.org/licenses/.

#' @rdname rFerns
#' @export
rFerns<-function(x,...)
 UseMethod("rFerns")

 #' @rdname rFerns
 #' @method rFerns formula
 #' @param formula alternatively, formula describing model to be analysed.
 #' @param data in which to interpret formula.
 #' @export
rFerns.formula<-function(formula,data=.GlobalEnv,...){
 #Convert formula into a data frame
 stats::terms.formula(formula,data=data)->t
 rx<-eval(attr(t,"variables"),data)
 apply(attr(t,"factors"),1,sum)>0->sel
 nam<-rownames(attr(t,"factors"))[sel]
 data.frame(rx[sel])->x; names(x)<-nam
 rx[[attr(t,"response")]]->y
 #Pass to the default method
 rFerns.default(x,y,...)
}

#' @rdname rFerns
#' @method rFerns matrix
#' @export
rFerns.matrix<-function(x,y,...){
 #If the input is matrix, data.frame it first
 rFerns.default(data.frame(x),y,...)
}


#' Classification with random ferns
#'
#' This function builds a random ferns model on the given training data.
#' @rdname rFerns
#' @method rFerns default
#' @param x Data frame containing attributes; must have unique names and contain only numeric, integer or (ordered) factor columns.
#' Factors must have less than 31 levels. No \code{NA} values are permitted.
#' @param y A decision vector. Must a factor of the same length as \code{nrow(X)} for ordinary many-label classification, or a logical matrix with each column corresponding to a class for multi-label classification.
#' @param depth The depth of the ferns; must be in 1--16 range. Note that time and memory requirements scale with \code{2^depth}.
#' @param ferns Number of ferns to be build.
#' @param importance Set to calculate attribute importance measure (VIM);
#' \code{"simple"} will calculate the default mean decrease of true class score (MDTS, something similar to Random Forest's MDA/MeanDecreaseAccuracy),
#' \code{"shadow"} will calculate MDTS and additionally MDTS of this attribute shadow, an implicit feature build by shuffling values within it, thus stripping it from information (which is slightly slower).
#' Shadow importance is useful as a reference to judge significance of a regular importance.
#' \code{"none"} turns importance calculation off, for a slightly faster execution.
#' For compatibility with pre-1.2 rFerns, \code{TRUE} will resolve to \code{"simple"} and \code{FALSE} to \code{"none"}.
#' Abbreviation can be used instead of a full value.
#' @param reportErrorEvery If set to a number larger than 0, current OOB error approximation will be printed every time a chunk of \code{reportErrorEvery} ferns is finished.
#' @param saveErrorPropagation Should the OOB error approximation be calculated after each ferns was created and saved?
#' Setting to \code{FALSE} may improve performance.
#' @param saveForest Should the model be saved? It must be \code{TRUE} if you want to use the model for prediction; however, if you are interested in importance or OOB error only, setting it to \code{FALSE} significantly improves memory requirements, especially for large \code{depth} and \code{ferns}.
#' @param consistentSeed PRNG seed used for shadow importance \emph{only}.
#' Must be either a 2-element integer vector or \code{NULL}, which corresponds to seeding from the default PRNG.
#' It should be set to the same value in two merged models to make shadow importance meaningful.
#' @param ... For formula and matrix methods, a place to state parameters to be passed to default method.
#' For the print method, arguments to be passed to \code{print}.
#' @return An object of class \code{rFerns}, which is a list with the  following components:
#' \item{model}{The built model; \code{NULL} if \code{saveForest} was \code{FALSE}.}
#' \item{oobErr}{OOB approximation of accuracy.
#' Ignores never-OOB-tested objects (see oobScores element).}
#' \item{importance}{The importance scores or \code{NULL} if \code{importance} was set to \code{"none"}.
#' In a first case it is a \code{data.frame} with two or three columns:
#' \code{MeanScoreLoss} which is a mean decrease of a score of a correct class when a certain attribute is permuted,
#' \code{Tries} which is number of ferns which utilised certain attribute, and, only when \code{importance} was set to \code{"shadow"},
#' \code{Shadow}, which is a mean decrease of accuracy for the correct class for a permuted copy of an attribute (useful as a baseline for normal importance).
#' The \code{rownames} are set and equal to the \code{names(x)}.}
#' \item{oobScores}{A matrix of OOB scores of each class for each object in training set.
#' Rows correspond to classes in the same order as in \code{levels(Y)}.
#' If the \code{ferns} is too small, some columns may contain \code{NA}s, what means that certain objects were never in test set.}
#' \item{oobPreds}{A vector of OOB predictions of class for each object in training set. Never-OOB-tested objects (see above) have predictions equal to \code{NA}.}
#' \item{oobConfusionMatrix}{Confusion matrix build from \code{oobPreds} and \code{y}.}
#' \item{timeTaken}{Time used to train the model (smaller than wall time because data preparation and model final touches are excluded; however it includes the time needed to compute importance, if it applies).
#' An object of \code{difftime} class.}
#' \item{parameters}{Numerical vector of three elements: \code{classes}, \code{depth} and \code{ferns}, containing respectively the number of classes in decision and copies of \code{depth} and \code{ferns} parameters.}
#' \item{classLabels}{Copy of \code{levels(Y)} after purging unused levels.}
#' \item{consistentSeed}{Consistent seed used; only present for \code{importance="shadow"}.
#' Can be used to seed a new model via \code{consistentSeed} argument.}
#' \item{isStruct}{Copy of the train set structure, required internally by predict method.}
#' @note The unused levels of the decision will be removed; on the other hand unused levels of categorical attributes will be preserved, so that they could be present in the data later predicted with the model.
#' The levels of ordered factors in training and predicted data must be identical.
#'
#' Do not use formula interface for a data with large number of attributes; the overhead from handling the formula may be significant.
#' @references Ozuysal M, Calonder M, Lepetit V & Fua P. (2010). \emph{Fast Keypoint Recognition using Random Ferns}, IEEE Transactions on Pattern Analysis and Machine Intelligence, 32(3), 448-461.
#'
#' Kursa MB (2014). \emph{rFerns: An Implementation of the Random Ferns Method for General-Purpose Machine Learning}, Journal of Statistical Software, 61(10), 1-13.
#' @author Miron B. Kursa
#' @examples
#' set.seed(77);
#' #Fetch Iris data
#' data(iris)
#' #Build model
#' rFerns(Species~.,data=iris)
#' ##Importance
#' rFerns(Species~.,data=iris,importance="shadow")->model
#' print(model$imp)
#' @export
#' @useDynLib rFerns random_ferns random_ferns_predict

rFerns.default<-function(x,y,depth=5,ferns=1000,importance="none",reportErrorEvery=0,saveErrorPropagation=FALSE,saveForest=TRUE,consistentSeed=NULL,...){
 #Stop on bad input
 depth<-as.integer(depth);
 ferns<-as.integer(ferns);
 stopifnot(length(depth)==1 && depth>0 && depth<=16)
 stopifnot(length(ferns)==1 && ferns>0)
 stopifnot(!any(is.na(y)))
 if(!is.data.frame(x)) stop("x must be a data frame.")
 if(is.na(names(x)) || any(duplicated(names(x)))) stop("Attribute names must be unique.")
 if(is.factor(y) && is.null(dim(y))){
  multi<-FALSE;
  if(length(y)!=nrow(x)) stop("Attributes' and decision's sizes must match.")
 }else{
  y<-as.matrix(y);
  if(is.logical(y) && length(dim(y))==2){
   multi<-TRUE;
   if(nrow(y)!=nrow(x)) stop("Attributes' and decision's sizes must match.")
  }else{
   stop("y must be a factor vector or a logical matrix.")
  }
 }
 if(!all(sapply(x,function(j) any(class(j)%in%c("numeric","integer","factor","ordered")))))
  stop("All attributes must be either numeric or factor.")
 if(any((sapply(x,function(a) ((length(levels(a))>30)&&(!is.ordered(a)))))->bad)){
  stop(sprintf(
   "Attribute(s) %s is/are unordered factor(s) with above 30 levels. Split or convert to ordered.",
   paste(names(x)[bad],collapse=", ")))
 }
 #Backward compatibility
 if(length(importance)!=1) stop("Wrong importance value.");
 if(identical(importance,FALSE)) importance<-"none";
 if(identical(importance,TRUE)) importance<-"simple";
 importance<-pmatch(importance,c("none","simple","shadow"));
 if(is.na(importance)) stop("Wrong importance value.");

 #Consistent seed setup
 if(importance==3){
  if(is.null(consistentSeed)){
   consistentSeed<-as.integer(sample(2^32-1,2,replace=TRUE)-2^31);
  }
  stopifnot(is.integer(consistentSeed));
  stopifnot(length(consistentSeed)==2);
 }else{
  if(!is.null(consistentSeed)){
   warning("Consistent seed is only useful with shadow importance; dropping.");
   consistentSeed<-NULL;
  }
 }

 if(multi && (importance>1)) stop("Importance is not yet supported for multi-label ferns.")

 if(reportErrorEvery<1 || reportErrorEvery>ferns) reportErrorEvery<-ferns+1
 saveOobErr<-ifelse(saveErrorPropagation,1,-1)*reportErrorEvery

 Sys.time()->before
 .Call(random_ferns,x,y,
  as.integer(depth[1]),
  as.integer(ferns[1]),
  as.integer(importance-1), #0->none, 1->msl, 2->msl+sha
  as.integer(saveOobErr),
  as.integer(saveForest),
  as.integer(multi),
  as.integer(consistentSeed))->ans
 after<-Sys.time()

 #Adjust C output with proper factor levels
 if(!is.null(ans$oobPreds)){
  ans$oobPreds<-factor(ans$oobPreds,
      levels=0:(length(levels(y))-1),
      labels=levels(y))
 }else{
  if(multi){
   ans$oobPreds<-t(ans$oobScores>0);
   colnames(ans$oobPreds)<-colnames(y);
  }
 }

 if(!multi){
  ans$classLabels<-levels(y)
 }else{
  ans$classLabels<-colnames(y)
 }
 if(saveForest){
  ans$isStruct<-list()
  lapply(x,levels)->ans$isStruct$predictorLevels
  sapply(x,is.integer)->ans$isStruct$integerPredictors
  sapply(x,is.ordered)->ans$isStruct$orderedFactorPredictors
 }

 if(!multi){
  table(Predicted=ans$oobPreds,True=y)->ans$oobConfusionMatrix
  if(is.null(ans$oobErr))
   ans$oobErr<-mean(ans$oobPreds!=y,na.rm=TRUE)
  ans$parameters<-c(classes=length(levels(y)),depth=depth,ferns=ferns);
  ans$type<-"class-many";
 }else{
  NULL->ans$oobConfusionMatrix
  if(is.null(ans$oobErr))
   ans$oobErr<-mean(rowSums(y!=ans$oobPreds))
  ans$oobPerClassError<-colMeans(ans$oobPreds!=y);
  ans$parameters<-c(classes=ncol(y),depth=depth,ferns=ferns);
  ans$type<-"class-multi";
 }

 if(!is.null(ans$importance)){
  if(importance==2){
   ans$importance<-data.frame(matrix(ans$importance,ncol=2))
   names(ans$importance)<-c("MeanScoreLoss","Tries")
  }
  if(importance==3){
   ans$importance<-data.frame(matrix(ans$importance,ncol=3))
   names(ans$importance)<-c("MeanScoreLoss","Shadow","Tries")
   ans$consistentSeed<-consistentSeed;
  }
  if(!is.null(names(x)))
   rownames(ans$importance)<-names(x)
 }

 #Calculate time taken by the calculation
 ans$timeTaken<-after-before

 class(ans)<-"rFerns"
 return(ans);
}

#' Prediction with random ferns model
#'
#' This function predicts classes of new objects with given \code{rFerns} object.
#' @method predict rFerns
#' @param  object Object of a class \code{rFerns}; a model that will be used for prediction.
#' @param x Data frame containing attributes; must have corresponding names to training set (although order is not important) and do not introduce new factor levels.
#' If this argument is not given, OOB predictions on the training set will be returned.
#' @param scores If \code{TRUE}, the result will contain score matrix instead of simple predictions.
#' @param ... Additional parameters.
#' @return Predictions.
#' If \code{scores} is \code{TRUE}, a factor vector (for many-class classification) or a logical data.frame (for multi-class classification) with predictions, else a data.frame with class' scores.
#' @author Miron B. Kursa
#' @examples
#' set.seed(77)
#' #Fetch Iris data
#' data(iris)
#' #Split into tRain and tEst set
#' iris[c(TRUE,FALSE),]->irisR
#' iris[c(FALSE,TRUE),]->irisE
#' #Build model
#' rFerns(Species~.,data=irisR)->model
#' print(model)
#'
#' #Test
#' predict(model,irisE)->p
#' print(table(
#'  Predictions=p,
#'  True=irisE[["Species"]]))
#' err<-mean(p!=irisE[["Species"]])
#' print(paste("Test error",err,sep=" "))
#'
#' #Show first OOB scores
#' head(predict(model,scores=TRUE))
#' @export
predict.rFerns<-function(object,x,scores=FALSE,...){
 #Validate input
 if(!("rFerns"%in%class(object))) stop("object must be of a rFerns class")
 if(is.null(object$model)&(!missing(x)))
  stop("This fern forest object does not contain the model.")
 scores<-as.logical(scores)[1];
 if(is.na(scores)) stop("Wrong value of scores; should be TRUE or FALSE.")

 iss<-object$isStruct
 if(is.null(iss)){
  #object is a v0.1 rFerns
  object$isStruct$predictorLevels<-object$predictorLevels
  object$isStruct$integerPredictors<-
   object$isStruct$orderedFactorPredictors<-
    rep(FALSE,length(iss$predictorLevels))
 }
 iss$predictorLevels->pL
 pN<-names(pL)
 multi<-identical(object$type,"class-multi");

 if(missing(x)){
  #There is no x; return the OOB predictions, nicely formatted
  if(scores){
   data.frame(t(object$oobScores))->ans;
   object$classLabels->names(ans)
   return(ans);
  }else{
   if(multi){
    return(data.frame(object$oobPreds));
   }else{
    return(object$oobPreds);
   }
  }
 }

 if(!identical(names(x),pN)){
  #Restore x state from training based on x's names
  if(!all(pN%in%names(x))){
   stop("Some training attributes missing in test.");
  }
  x[,pN]->x
 }

 #Fail for NAs in input
 if(any(is.na(x))) stop("NAs in predictors.")

 for(e in 1:ncol(x))
  if(is.null(pL[[e]])){
   if(iss$integerPredictors[e]){
    if(!("integer"%in%class(x[,e]))) stop(sprintf("Attribute %s should be integer.",pN[e]))
   }else{
    if(!("numeric"%in%class(x[,e]))) stop(sprintf("Attribute %s should be numeric.",pN[e]))
   }
  }else{
   if(iss$orderedFactorPredictors[e]){
    #Check if given attribute is also ordered
    if(!is.ordered(x[,e])) stop(sprintf("Attribute %s should be an ordered factor.",pN[e]))
    #Convert levels
    if(!identical(levels(x[,e]),pL[[e]]))
     stop(sprintf("Levels of %s does not match those from training (%s).",pN[e],paste(pL[[e]],collapse=", ")))
   }else{
    #Convert factor levels to be compatible with training
    if(!identical(levels(x[,e]),pL[[e]]))
     x[,e]<-factor(x[,e],levels=pL[[e]])
    #In case of mismatch, NAs will appear -- catch 'em and fail
    if(any(is.na(x[,e]))) stop(sprintf("Levels of %s does not match those from training (%s).",pN[e],paste(pL[[e]],collapse=", ")))
   }
  }

 #Prediction itself
 Sys.time()->before
 .Call(random_ferns_predict,x,
  object$model,
  as.integer(object$parameters["depth"]),
  as.integer(object$parameters["ferns"]),
  as.integer(length(object$classLabels)),
  as.integer(scores),as.integer(multi))->ans
 after<-Sys.time()

 if(scores){
  ans<-data.frame(matrix(ans,ncol=length(object$classLabels),byrow=TRUE)/
   object$parameters["ferns"])
  object$classLabels->names(ans)
 }else{
  if(!multi){
   ans<-factor(ans,levels=0:(length(object$classLabels)-1),
    labels=object$classLabels)
  }else{
   ans<-data.frame(matrix(ans,ncol=length(object$classLabels),byrow=TRUE)>0)
   object$classLabels->names(ans)
  }
 }

 #Store the timing
 attr(ans,"timeTaken")<-after-before

 return(ans)
}

#' @method print rFerns
#' @export
print.rFerns<-function(x,...){
 #Pretty-print rFerns output
 cat(sprintf("\n Forest of %s %sferns of a depth %s.\n\n",
  x$parameters["ferns"],
  ifelse(identical(x$type,"class-multi"),"multi-class ",""),
  x$parameters["depth"]))
 if(identical(x$type,"class-multi")){
  if(!is.null(x$oobErr))
   cat(sprintf(" OOB Hamming distance %0.3f for %d classes.\n",
    utils::tail(x$oobErr,1),
    x$parameters["classes"]))
   if(!is.null(x$oobPerClassError)){
    cat(" Per-class error rates:\n");
    print(x$oobPerClassError)
   }
 }else{
  if(!is.null(x$oobErr))
   cat(sprintf(" OOB error %0.2f%%;",utils::tail(x$oobErr,1)*100))
  if(!is.null(x$oobConfusionMatrix)){
   cat(" OOB confusion matrix:\n");
   print(x$oobConfusionMatrix)
  }
 }
 if(!is.null(x$oobScores) && any(is.na(x$oobScores)))
  cat(" Note: forest too small to provide good OOB approx.\n");

 return(invisible(x));
}
