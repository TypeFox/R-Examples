#' Merge two random ferns models
#'
#' This function combines two compatible (same decision, same training data structure and same depth) models into a single ensemble.
#' It can be used to distribute model training, perform it on batches of data, save checkouts or precisely investigate its course.
#' @method merge rFerns
#' @param x Object of a class \code{rFerns}; a first model to be merged.
#' @param y Object of a class \code{rFerns}; a second model to be merged.
#' Can also be \code{NULL}, \code{x} is immediately returned in that case.
#' Has to have be built on the same kind of training data as \code{x}, with the same depth.
#' @param dropModel If \code{TRUE}, model structure will be dropped to save size.
#' This disallows prediction using the merged model, but retains importance and OOB approximations.
#' @param ignoreObjectConsistency If \code{TRUE}, merge will be done even if both models were built on a different sets of objects.
#' This drops OOB approximations.
#' @param trueY Copy of the training decision, used to re-construct OOB error and confusion matrix.
#' Can be omitted, OOB error and confusion matrix will disappear in that case; ignored when \code{ignoreObjectConsistency} is \code{TRUE}.
#' @param ... Ignored, for S3 gerneric/method consistency.
#' @return An object of class \code{rFerns}, which is a list with the  following components:
#' \item{model}{The merged model in case both \code{x} and \code{y} had model structures included and \code{dropModel} was \code{FALSE}.
#' Otherwise \code{NULL}.}
#' \item{oobErr}{OOB approximation of accuracy, if can be computed.
#' Namely, when \code{oobScores} could be and \code{trueY} is provided.}
#' \item{importance}{The merged importance scores in case both \code{x} and \code{y} had importance calculated.
#' Shadow importance appears only if both models had it enabled.}
#' \item{oobScores}{OOB scores, if can be computed; namely if both models had it calculated and \code{ignoreObjectConsistency} was not used.}
#' \item{oobPreds}{A vector of OOB predictions of class for each object in training set, if can be computed.}
#' \item{oobConfusionMatrix}{OOB confusion matrix, if can be computed.
#' Namely, when \code{oobScores} could be and \code{trueY} is provided.}
#' \item{timeTaken}{Time used to train the model, calculated as a sum of training times of \code{x} and \code{y}.}
#' \item{parameters}{Numerical vector of three elements: \code{classes}, \code{depth} and \code{ferns}.}
#' \item{classLabels}{Copy of \code{levels(Y)} after purging unused levels.}
#' \item{isStruct}{Copy of the train set structure.}
#' \item{merged}{Set to \code{TRUE} to mark that merging was done.}
#' @note In case of different training object sets were used to build the merged models, merged importance is calculated but mileage may vary; for substantially different sets it may become biased.
#' Your have been warned.
#'
#' Shadow importance is only merged when both models have shadow importance and the same \code{consistentSeed} value; otherwise shadow importance would be biased down.
#'
#' The order of objects in \code{x} and \code{y} is not important; the only exception is merging with \code{NULL}, in which case \code{x} must be an \code{rFerns} object for R to use proper merge method.
#' @author Miron B. Kursa
#' @examples
#' set.seed(77)
#' #Fetch Iris data
#' data(iris)
#' #Build models
#' rFerns(Species~.,data=iris)->modelA
#' rFerns(Species~.,data=iris)->modelB
#' modelAB<-merge(modelA,modelB);
#' print(modelA);
#' print(modelAB);
#' @export
merge.rFerns<-function(x,y,dropModel=FALSE,ignoreObjectConsistency=FALSE,trueY=NULL,...){
 stopifnot(inherits(x,'rFerns')); #Tautology thanks to object dispatch
 if(is.null(y)) return(x);
 stopifnot(inherits(y,'rFerns'))
 stopifnot(identical(x$isStruct,y$isStruct));
 if(!ignoreObjectConsistency){
  stopifnot(identical(dim(x$oobScores),dim(y$oobScores)));
 }else{
  x$oobScores<-NULL;
  y$oobScores<-NULL;
  trueY<-NULL;
 }
 stopifnot(identical(x$classLabels,y$classLabels));
 stopifnot(identical(x$type,y$type));
 stopifnot(identical(x$parameters[-3],y$parameters[-3]))

 #Initiate core structure
 ans<-list(
  isStruct=x$isStruct,
  type=x$type,
  classLabels=x$classLabels,
  merged=TRUE);

 #Merge model fields
 if(is.null(x$model)||is.null(y$model)||dropModel){
  ans$model<-NULL;
 }else{
  ans$model<-list(
   splitAttIdxs=c(x$model$splitAttIdxs,y$model$splitAttIdxs),
   threReal=c(x$model$threReal,y$model$threReal),
   threInteger=c(x$model$threInteger,y$model$threInteger),
   scores=c(x$model$scores,y$model$scores)
  );
 }

 #OOB scores is summed with #ferns-derieved weights, then OOB elements re-generated
 #When nrow x!=nrow y, we assume it means different batches were used and oobScores make no sense anymore
 if(!is.null(x$oobScores)&&!is.null(y$oobScores)){
  ans$oobScores<-x$oobScores+y$oobScores;
  ans$oobPreds<-factor(x$classLabels)[apply(ans$oobScores,2,
   function(x){
    which.max(x)->l;
    if(length(l)!=1) return(NA);
    return(l);
   }
   )];
  if(!is.null(trueY)){
    #TODO: Multilabel!
   stopifnot(is.factor(trueY));
   stopifnot(identical(levels(trueY),ans$classLabels));
   #OOB error propagation is always squashed into final OOB error
   ans$oobErr<-mean(trueY!=ans$oobPreds);
   ans$oobConfusionMatrix<-table(Predicted=ans$oobPreds,True=trueY)
  }
 }

 #Importance
 if(!is.null(x$importance)&&!is.null(y$importance)){
  mergeWeight<-c(x$parameters["ferns"],y$parameters["ferns"]);
  mergeWeight<-mergeWeight/sum(mergeWeight);
  ans$importance<-x$importance;
  if(!is.null(ans$importance$MeanScoreLoss)&&!is.null(y$importance$MeanScoreLoss)){
   ans$importance$MeanScoreLoss<-
    mergeWeight[1]*x$importance$MeanScoreLoss+
    mergeWeight[2]*y$importance$MeanScoreLoss;
  }else{
   ans$importance$MeanScoreLoss<-NULL;
  }
  if(!is.null(ans$importance$Tries)&&!is.null(y$importance$Tries)){
   ans$importance$Tries<-x$importance$Tries+y$importance$Tries;
  }else{
   ans$importance$Tries<-NULL;
  }
  if(!is.null(ans$importance$Shadow)&&!is.null(y$importance$Shadow)&&identical(x$consistentSeed,y$consistentSeed)){
   ans$importance$Shadow<-
    mergeWeight[1]*x$importance$Shadow+
    mergeWeight[2]*y$importance$Shadow;
  }else{
   ans$importance$Shadow<-NULL;
  }
  if(!is.null(ans$importance$Hits)&&!is.null(y$importance$Hits)&&identical(x$consistentSeed,y$consistentSeed)){
   ans$importance$Hits<-x$importance$Hits+y$importance$Hits;
  }else{
   ans$importance$Hits<-NULL;
  }
 }

 #Parameters
 ans$parameters<-c(
  x$parameters["classes"],
  x$parameters["depth"],
  x$parameters["ferns"]+y$parameters["ferns"]);

 #Time taken
 ans$timeTaken<-x$timeTaken+y$timeTaken;

 class(ans)<-"rFerns";
 ans
}
