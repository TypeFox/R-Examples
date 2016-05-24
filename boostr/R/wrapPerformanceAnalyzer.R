#' @family Wrapper Generators
#' 
#' @title Create a boostr compatible wrapper for a performance analyzer.
#' 
#' @description Use provided metadata on a given performance analyzer to create
#' a boostr compatible wrapper.
#' 
#' @param analyzePerformance a function to analyze the performance of an estimator
#' @param analyzerInputPreds a string indicating the name of the argument in
#' \code{analyzePerformance}'s signature that represents the estimator's
#' predictions.
#' @param analyzerInputResponse a string indicating the name of the argument in
#' \code{analyzePerformance}'s signature that represents the true response 
#' associated with the estimator's predictions.
#' @param analyzerInputOOBObs a string indiciating the name of the argument in
#' \code{analyzePerformance}'s signature that represents the vector of indices
#' indicating which observations were out-of-bag.
#' @param .verbose a boolean indicating if warnings should be displayed or not.
#' 
#' @details
#' Since "performance" is a subjective thing, the requirements for a function to
#' be wrappable by \code{\link{wrapPerformanceAnalyzer}} are that they accept
#' predictions, true responses, and a vector of indices for out-of-bag
#' observations. After each iteration of the ensemble building phase in 
#' \code{\link{boostBackend}}, these three objects are fed to a performance 
#' analyzer. The output of the performance analyze is stored in the 
#' \code{estimatorPerformance} attribute of the object returned by 
#' \code{\link{boostBackend}}. 
#' 
#' @template performanceAnalyzers
#' 
#' @return A function (wrapper around \code{analyzePerformance}) which is also
#' a '\code{performanceAnalyzer}' object. The function's signature is
#' \code{(prediction, response, oobObs, ...)} and it's output preserves the
#'  output of \code{analyzePerformance}. Hence, the wrapper is a boostr
#' compatible performance analyzer. 
wrapPerformanceAnalyzer <- function(analyzePerformance, 
                                    analyzerInputPreds = "prediction",
                                    analyzerInputResponse = "response",
                                    analyzerInputOOBObs = "oobObs",
                                    .verbose = FALSE) {
  func <- addDots(analyzePerformance, .verbose)
  wrapper <- function(prediction, response, oobObs, ...) {
    
    funcArgs <- list(prediction = prediction,
                     response = response,
                     oobObs = oobObs)
    
    names(funcArgs) <- c(analyzerInputPreds,
                         analyzerInputResponse,
                         analyzerInputOOBObs)
    
    do.call(func, c(funcArgs, list(...)))
  }
  
  class(wrapper) <- c("performanceAnalyzer", class(wrapper))
  return(wrapper)
}