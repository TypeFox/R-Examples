#' @rdname boostWith
#' 
#' @title Boostr implemented versions of arc-fs, arc-x4 and AdaBoost.
#' @description
#' Perform the Boost algorithm for the algorithms arc-fs, arc-x4, and AdaBoost.
#' 
#' @param x a list with entries '\code{train}' and '\code{predict}' or a 
#' function that satisfies the definition of an estimation procedure given 
#' below. The list input will invoke a call to 
#' \code{\link{buildEstimationProcedure}}. Function input will invoke a call to
#' \code{\link{wrapProcedure}}, unless the function inherits from 
#' '\code{estimationProcedure}'. In either event, metadata may be required to
#' properly wrap \code{x}. See the appropriate help documentation.
#' @param B number of iterations of boost to perform.
#' @param data a data.frame of matrix to act as the learning set. The columns
#'  are assumed to be ordered such that the response variable in the first
#'  column and the remaining columns as the predictors. As a convenience,
#'  \code{\link{boostBackend}} comes with a switch, \code{.formatData} 
#'  (defaulted to \code{TRUE}) which will look for  an argument named 
#'  \code{formula} inside \code{.procArgs} and use the value of
#'  \code{formula} to format \code{data}. If you don't want this to happen, 
#'  or if the data is already properly formatted, include 
#'  \code{.formatData=FALSE} in \code{metadata}.
#' @param .procArgs a named list of arguments to pass to the estimation 
#' procedure.
#' If \code{x} is a list, \code{.procArgs} is a named list of lists with
#' entries \code{.trainArgs} and \code{.predictArgs} and each list is a 
#' named list of arguments to pass to \code{x$train} and \code{x$predict},
#' respectively. If \code{x} is a function, \code{.procArgs} is a named list 
#' of arguments to pass to \code{x}, in addition to \code{data} and 
#' \code{weights}. See 'Examples' below.
#' @param initialWeights a vector of weights used for the first iteration of
#' the ensemble building phase of Boost.
#' @param analyzePerformance a function which accepts an estimator's 
#' predictions and the true responses to said predictions (among other
#' arguments) and returns a list of values. If no function is provided, 
#' \code{\link{defaultOOBPerformanceAnalysis}} is used.
#'  See \code{\link{wrapPerformanceAnalyzer}} for metadata that may
#' need to be passed to make \code{analyzePerformance} compatible with the
#' boostr framework.
#' @param metadata a named list of additional arguments to be passed to 
#' \code{\link{wrapProcedure}}, \code{\link{buildEstimationProcedure}},
#' \code{\link{wrapPerformanceAnalyzer}} and/or \code{\link{boostBackend}}. 
#' @param .boostBackendArgs a named list of additional arguments to pass to
#' \code{\link{boostBackend}}.
#' 
#' @details
#' These functions call \code{\link{boost}} with the appropriate reweighters,
#' aggregators, and metadata.
#' 
#' @return a "\code{boostr}" object that is the output of 
#' \code{\link{boostBackend}}. 
#' @export
boostWithArcFs <-
  function(x, B, data, .procArgs=NULL, metadata=NULL,
           initialWeights=rep.int(1, nrow(data))/nrow(data),
           analyzePerformance = defaultOOBPerformanceAnalysis,
           .boostBackendArgs=NULL) {
    
  boost(x=x, B=B, data=data,.procArgs=.procArgs,
        metadata=metadata,
        initialWeights=initialWeights,
        analyzePerformance=analyzePerformance,
        reweighter=arcfsReweighter,
        aggregator=arcfsAggregator,
        .boostBackendArgs=.boostBackendArgs)
}

#' @rdname boostWith
#' @export
boostWithArcX4 <-
  function(x, B, data, .procArgs=NULL, metadata=NULL, 
           initialWeights=rep.int(1, nrow(data))/nrow(data),
           analyzePerformance = defaultOOBPerformanceAnalysis,
           .boostBackendArgs = NULL) {
    
    .boostBackendArgs <- c(.boostBackendArgs, list(.reweighterArgs=list(m=0)))
    
    boost(x=x, B=B, data=data,.procArgs=.procArgs,
          metadata=metadata,
          initialWeights=initialWeights,
          analyzePerformance=analyzePerformance,
          reweighter=arcx4Reweighter,
          aggregator=arcx4Aggregator,
          .boostBackendArgs=.boostBackendArgs)
  }

#' @rdname boostWith
#' @export
boostWithAdaBoost <-
  function(x, B, data, .procArgs=NULL, metadata=NULL,
           initialWeights=rep.int(1, nrow(data))/nrow(data),
           analyzePerformance = defaultOOBPerformanceAnalysis,
           .boostBackendArgs = NULL) {
    
    boost(x=x, B=B, data=data,.procArgs=.procArgs,
          metadata=metadata,
          initialWeights=initialWeights,
          analyzePerformance=analyzePerformance,
          reweighter=adaboostReweighter,
          aggregator=adaboostAggregator,
          .boostBackendArgs=.boostBackendArgs)
  }