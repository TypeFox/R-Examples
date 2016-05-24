#' @family Wrapper Generators
#' 
#' @title Create a boostr compatible wrapper for an estimation procedure.
#' 
#' @description Use provided metadata on a given estimation procedure to create
#' a boostr compatible wrapper. See section below for more details on estimation
#' procedures.
#'
#' @param proc a function that obeys the definition of an estimation procedure
#' as defined in the white paper. Generally, \code{proc} must be a function
#' which learns some model and consequently returns an estimator that uses the
#' learned model. See below.
#' @param learningSet  a string indicating the name of the argument in
#' \code{proc}'s signature that passes the data to be used inside \code{proc}. 
#' @param predictionSet a string indicating the name of the argument in
#' \code{predict}'s signature that indicates the observation to predicate
#' responses for.
#' 
#' @template estimationProcedures
#' 
#' @return An '\code{estimationProcedure}' object whose signature and whose
#' output's signature are compatible with boostr. Explicitly, the arguments of
#' the wrapper are
#' \item{data}{the data that \code{proc} will use to build a model.}
#' \item{...}{any additional arguments necessary for \code{proc} to make its
#' model.}
#' and the returned closure from the wrapper has arguments
#' \item{newdata}{the data that \code{proc}'s output will predict responses for.}
#' \item{.estimatorArgs}{a named list of any additional arguments that need to
#' be passed to \code{proc}'s output.}
#' @examples
#' \dontrun{
#' library(class)
#' 
#' # an estimation procedure outputs a estimator
#' 
#' knnProc <- function(formula, traindata, k) {
#'   df <- model.frame(formula=formula, data=traindata)
#'   function(testdata, prob=FALSE) {
#'     knn(train=df[, -1], test=testdata, cl=df[, 1], prob=prob, k=k) 
#'   }
#' }
#' 
#' boostrKNN <- wrapProcedure(knnProc,
#'                           learningSet="traindata",
#'                           predictionSet="testdata")
#' }

wrapProcedure <- function(proc, learningSet="data", predictionSet="newdata") {
  # build function factory (Phi)
  f <- function(data, ...) {
    # build estimator according to proc and metadata
    .procArgs <- c(list(...), list(data=data))
    names(.procArgs)[length(.procArgs)] <- learningSet
    
    estimator <- do.call(proc, .procArgs)
    
    # return a wrapper of estimator with a reformatted signature.
    function(newdata, .estimatorArgs=NULL) {
      .estimatorArgs <- c(list(newdata=newdata), .estimatorArgs)
      names(.estimatorArgs)[1] <- predictionSet
      
      do.call(estimator, .estimatorArgs)
    }
  }
  
  class(f) <- c("estimationProcedure", class(f))
  f
}
