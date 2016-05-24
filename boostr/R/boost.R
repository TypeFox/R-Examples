# helper functions for below.
extractReweighterMetadata <- function(metadata) {
  metadata[names(metadata) %in% names(formals(wrapReweighter))]
}

extractAggregatorMetadata <- function(metadata) {
  metadata[names(metadata) %in% names(formals(wrapAggregator))]
}

extractAnalyzerMetadata <- function(metadata) {
  metadata[names(metadata) %in% names(formals(wrapPerformanceAnalyzer))]
}

extractBoostBackendArgs <- function(metadata) {
  metadata[names(metadata) %in% names(formals(boostBackend))]
}

# a function factory that abstracts the procedures in
# boost.list and boost.function
boostAccordingToWrapperBuilder <- function(wrapperBuilder) {
  
  function(x, B, reweighter, aggregator, data, .procArgs, metadata,
           initialWeights, analyzePerformance, .boostBackendArgs) {
    
    extractProcedureMetadata <- function(metadata) {
        metadata[names(metadata) %in% names(formals(wrapperBuilder))] 
    }
    
    # helper function intended solely for its side-effects
    wrapReweighterAggregatorAndAnalyzer <- function(.reweighter, .aggregator,
                                                 .analyzer, metadata) {
      
      if (!inherits(.reweighter, 'reweighter')) {
        .reweighterMetadata <- extractReweighterMetadata(metadata)
        reweighter <<- do.call(wrapReweighter,
                               c(reweighter=.reweighter, .reweighterMetadata))
      }
      
      if (!inherits(.aggregator, 'aggregator')) {
        .aggregatorMetadata <- extractAggregatorMetadata(metadata)
        aggregator <<- do.call(wrapAggregator, 
                               c(aggregator=.aggregator,.aggregatorMetadata))
      }
      
      if (!is.null(.analyzer) && !inherits(.analyzer, "performanceAnalyzer")) {
        .analyzerMetadata <- extractAnalyzerMetadata(metadata)
        analyzePerformance <<- do.call(wrapPerformanceAnalyzer, 
                               c(analyzePerformance=.analyzer,
                                 .analyzerMetadata))
      }
    }
    
    .procMetadata <- extractProcedureMetadata(metadata)
    
    estimationProc <- do.call(wrapperBuilder, c(x, .procMetadata))
    
    wrapReweighterAggregatorAndAnalyzer(reweighter, aggregator,
                                     analyzePerformance, metadata)
    
    .boostBackendArgs <- 
      c(B=B, reweighter=reweighter, aggregator=aggregator, proc=estimationProc,
        list(.procArgs=.procArgs, data=data, initialWeights=initialWeights),
        .boostBackendArgs, analyzePerformance=analyzePerformance)
    
    out <- do.call(boostBackend, .boostBackendArgs)
    if (!inherits(out, "boostr")) class(out) <- c('boostr', class(out))
    return(out)
  }
}

#' @title Boost an Estimation Procedure with a Reweighter and an Aggregator.
#' @description Boost an \emph{estimation procedure} and analyze individual
#' estimator performance using a \emph{reweighter}, \emph{aggregator}, and 
#' some \emph{performance analyzer}. 
#' 
#' @param B number of iterations of boost to perform.
#' @param x a list with entries '\code{train}' and '\code{predict}' or a 
#' function that satisfies the definition of an estimation procedure given 
#' below. The list input will invoke a call to 
#' \code{\link{buildEstimationProcedure}}. Function input will invoke a call to
#' \code{\link{wrapProcedure}}, unless the function inherits from 
#' '\code{estimationProcedure}'. In either event, metadata may be required to
#' properly wrap \code{x}. See the appropriate help documentation.
#' @param reweighter A reweighter, as defined below. If the function does not
#' inherit from '\code{reweighter}', a call to \code{\link{wrapReweighter}}
#' will be made. See \code{\link{wrapReweighter}} to determine what metadata,
#' if any, you may need to pass for the wrapper to be \code{boostr} compatible
#' @param aggregator An aggregator, as defined below. If the function does not
#'  inherit from '\code{aggregator}' a call to \code{\link{wrapAggregator}}
#'  will be made to build a boostr compatible wrapper. See 
#'  \code{\link{wrapAggregator}} to determine if any  metadata needs to be
#'  passed in for this to be successful.
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
#' @param metadata a named list of arguments to be passed to 
#' \code{\link{wrapProcedure}}, \code{\link{buildEstimationProcedure}},
#'  \code{\link{wrapReweighter}}, \code{\link{wrapAggregator}}, and/or
#'  \code{\link{wrapPerformanceAnalyzer}}. 
#' @param .boostBackendArgs a named list of additional arguments to pass to
#' \code{\link{boostBackend}}. 
#'  
#' @details This function is a designed to be an interface between the user
#' and \code{\link{boostBackend}} when \code{x}, \code{reweighter},
#' \code{aggregator} and/or \code{analyzePerformance} are valid input to 
#' the Boost algorithm, but do not have \code{boostr} compatible signatures.
#' Hence, \code{boost} calls the appropriate wrapper function (with the
#' relevant information from \code{metadata}) to convert user supplied
#' functions into \code{boostr} compatible functions.
#' 
#' @rdname boost
#' 
#' @family boostr wrappers
#' @family performance analyzers
#' @family aggregators
#' @family reweighters
#' 
#'  @export
#'  @return
#'  a '\code{boostr}' object which is returned from \code{\link{boostBackend}}.
#'  This object is a function of a single input
#'  \item{newdata}{a data.frame or matrix whose columns should probably be in
#'   the same order as the columns of the data each of the constituent
#'   estimators was trained on.}
#'  The return value of this function is a prediction for each row in
#'  \code{newdata}.
#'  
#'  See \code{\link{boostBackend}} for more details on "\code{boostr}" objects.
#'  
#' @example inst/examples/test-boost_examples.R
boost <- function(x, B, reweighter, aggregator, data,
                  .procArgs=NULL, metadata=NULL,
                  initialWeights=rep.int(1, nrow(data))/nrow(data),
                  analyzePerformance = defaultOOBPerformanceAnalysis,
                  .boostBackendArgs = NULL) {
  UseMethod("boost")
}

#' @rdname boost
#' @export
boost.list <- function(x, B, reweighter, aggregator, data,
                       .procArgs=NULL, metadata=NULL,
                       initialWeights=rep.int(1, nrow(data))/nrow(data),
                       analyzePerformance = defaultOOBPerformanceAnalysis,
                       .boostBackendArgs = NULL) {
  f <- boostAccordingToWrapperBuilder(buildEstimationProcedure)
  f(x=x, B=B, reweighter=reweighter, aggregator=aggregator,
    data=data, .procArgs=.procArgs, metadata=metadata,
    initialWeights=initialWeights, .boostBackendArgs = .boostBackendArgs,
    analyzePerformance=analyzePerformance)
}


#' @rdname boost
#' @export
boost.function <- function(x, B, reweighter, aggregator, data,
                           .procArgs=NULL, metadata=NULL,
                           initialWeights=rep.int(1, nrow(data))/nrow(data),
                           analyzePerformance= defaultOOBPerformanceAnalysis,
                           .boostBackendArgs = NULL) {
  
  f <- boostAccordingToWrapperBuilder(wrapProcedure)
  
  f(x=x, B=B, reweighter=reweighter, aggregator=aggregator, data=data,
    .procArgs=.procArgs, metadata=metadata, initialWeights=initialWeights,
    analyzePerformance=analyzePerformance, .boostBackendArgs = .boostBackendArgs)
}