#' @include 05PredictionControlClass.R
NULL

# setOldClass("C5.0")

#' PreProCombClass
#'
#' PreproCombClass is a container for output of executing preprocomb function.
#' It is an interface for extension of the system.
#' @slot rawall (data frame) all results
#' @slot catclassification (data frame) classification accuracy categorized, "high" is more than 80 percent quantile value
#' @slot allclassification (data frame) classification accuracy, means and standard deviations
#' @slot bestclassification (data frame) best classification accuracy combinations
#' @slot allclustering (data frame) hopkins statistics values
#' @slot bestclustering (data frame) best hopkins statistics combinations
#' @slot alloutliers (data frame) ORH outlier score for 95 percent quantile value
#' @details Extensions can included analysis of catclassification by association rules or modelling
#' @export

setClass("PreProCombClass", representation(rawall="data.frame", catclassification="data.frame", allclassification="data.frame", bestclassification="data.frame", allclustering="data.frame", bestclustering="data.frame", alloutliers="data.frame"))

#' preprocomb
#'
#' preprocomb the main function of programmatic use. It executes the computation of
#' classification accuracy, hopkins statistic and ORH outlier score. See also
#' PreProCombClass.
#
#' @param models (character) vector of models (names of models as defined in package caret)
#' @param gridclassobject (GridClass) object representing the grid of combinations
#' @param nholdout (integer) number of holdout rounds for predictive classification, must be two or more, defaults to two
#' @param searchmethod (character) defaults to "exhaustive" full blind search, "random" search 20 percent of grid, "grid" grid search 10 percent
#' @param predict (boolean) compute predictions, defaults to TRUE
#' @param cluster (boolean) compute clustering tendency, defaults to FALSE
#' @param outlier (boolean) compute outlier tendency, defaults to FALSE
#' @return a PreProCombClass object
#' @details caret messages will be displayed during processing
#' @examples
#' ## modifiediris <- droplevels(iris[-c(1:60),])
#' ## grid <- setgrid(phases=c("outliers", "scaling"), data=modifiediris)
#' ## library(kernlab)
#' ## result <- preprocomb(models=c("svmRadial"), grid=grid, nholdout=2, search="exhaustive")
#' ## result@@allclassification
#' ## result@@allclustering
#' ## result@@alloutliers
#' ## result@@rawall
#' ## result@@catclassification
#' ##
#' ## newphases <- c("outlier", "smoothing", "scaling", "selection", "sampling")
#' ## newmodels <- c("knn", "rf", "svmRadial")
#' ## grid1 <- setgrid(phases=newphases, data=modifiediris)
#' ## result1 <- preprocomb(models=newmodels, grid=grid1, nholdout=2, search="grid")
#' @export

preprocomb <- function(models="knn", gridclassobject, nholdout=2, searchmethod="exhaustive", predict=TRUE, cluster=FALSE, outlier=FALSE){

  predictors <- models

  # prediction control: TO BE ELIMINATED
  if(class(predictors)!="character"){stop("The argument predictors must a character vector.")}

  if(class(gridclassobject)!="GridClass"){stop("The argument grid must be a GridClass object.")}

  supportedsearches <- c("exhaustive", "random", "grid")
  if (!searchmethod %in% supportedsearches) {stop("Argument 'search' must on one of the following: 'exhaustive', 'random', 'grid")}

  if (class(predict)!="logical"){stop("Argument 'predict' must a logical")}
  if (class(cluster)!="logical"){stop("Argument 'cluster' must a logical")}
  if (class(outlier)!="logical"){stop("Argument 'outlier' must a logical")}

  out <- combinationevaluation(predictors, gridclassobject, nholdout, searchmethod, predict, cluster, outlier)

  preprocombclassobject <- new("PreProCombClass")

  # all classification accuracies

  preproout <- data.frame(matrix(paste(round(as.matrix(out[[1]]),2), round(as.matrix(out[[2]]),2), sep="+-"), nrow=nrow(out[[1]])))
  colnames(preproout) <- c(predictors, "ALL_MEAN")
  preprocombclassobject@allclassification <- data.frame(as.data.frame(out[[5]]), preproout)

  # best best classification accuracies

  bestaccuracies <-  head(preprocombclassobject@allclassification[order(preprocombclassobject@allclassification$ALL_MEAN, decreasing=TRUE),])
  preprocombclassobject@bestclassification <- bestaccuracies

  # raw data

  rawall <- data.frame(out[[5]], out[[1]], out[[2]], out[[3]], out[[4]])
  colnames(rawall)[(ncol(out[[5]])+1):ncol(rawall)] <- c(paste(c(predictors, "ALL_MEAN"), "Mean", sep=""), paste(c(predictors, "ALL_MEAN"), "SD", sep=""), "Hopkins", "Orh_skewness")
  preprocombclassobject@rawall <- rawall

  # raw categorical

  if (predict==TRUE){
  tempout <- data.frame(out[[5]], out[[1]][ncol(out[[1]])])
  colnames(tempout)[ncol(tempout)] <- "ALL_MEAN"
  cutpoint <- quantile(tempout$ALL_MEAN, .80)
  tempout$target <- cut(tempout$ALL_MEAN, breaks=c(-Inf, cutpoint, Inf), labels=c("low", "high"))
  tempout <- tempout[,-(ncol(tempout)-1)]
  preprocombclassobject@catclassification <- tempout
  }

  # by clustering tendency

  preprocombclassobject@allclustering <- data.frame(out[[5]], hopkins=round(out[[3]],2))

  bestclustering <- head(preprocombclassobject@allclustering[order(preprocombclassobject@allclustering$hopkins, decreasing=TRUE),])
  preprocombclassobject@bestclustering <- bestclustering

  # by outliers

  preprocombclassobject@alloutliers <- data.frame(out[[5]], Orh_skewness=round(out[[4]],2))

  return(preprocombclassobject)
}

