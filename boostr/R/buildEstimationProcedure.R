#' @family Wrapper Generators
#' @family procedure wrappers
#' 
#' @title Build a boostr compatible estimation procedure. 
#' 
#' @description A convenience function which builds a boostr compatible estimation
#' procedure from functions \code{train} and \code{predict}. 
#'
#' @param train a function that learns from data to produce a model
#' @param predict a function that leverages the model from \code{train} to 
#' generate predictions from new data.
#' @param learningSet a string indicating the name of the argument in
#' \code{train}'s signature that passes data inside \code{train}. 
#' @param predictionSet a string indicating the name of the argument in
#' \code{predict}'s signature that indicates the observation to predicate
#' responses for.
#' @param modelName a string indicating the name of the argument in
#' \code{predict}'s signature that passes the model from \code{train} inside
#' \code{predict}. 
#' 
#' @template estimationProcedures
#' 
#' @section Warning:
#'  This function makes the fundamental assumption that the design-pattern
#' linking \code{train} and {predict} is the common \code{train}-\code{predict}
#' pattern found in the design of \code{SVM} in the examples. If this is not the
#' case, you'll want to build assemble your procedure manually and call
#' \code{\link{wrapProcedure}} instead.
#' 
#' 
#' @return An '\code{estimationProcedure}' object which is compatible with the
#' boostr framework. Meaning, the output is a function factory  which accepts
#' arguments
#' \item{data}{the data to be passed to \code{train}.}
#' \item{.trainArgs}{a list of arguments to be passed to \code{train}, in
#'  addition to \code{data}. If the order of arguments in \code{train} is
#'  important, you'll need to respect that order inside \code{.trainArgs}.}
#' \item{.predictArgs}{a list of arguments to pass to \code{predict} in 
#' addition to \code{modelName} and \code{predictionSet}. If the order
#' of these arguments matters, respect that order in \code{.predictArgs}.}
#' 
#' and returns a closure with arguments
#' \item{newdata}{the data whose response variable is to be estimated.}
#' \item{.predictArgs}{a list of arguments to pass to \code{predict} in 
#' addition to \code{modelName} and \code{predictionSet}. This is defaulted
#' to the value of \code{.predictArgs} passed to the parent function, however
#' access to this has been granted as a convenience to the user. Again, if the
#' order of these arguments matters, respect that order in \code{.predictArgs}.}
#' 
#'  
#' @examples
#' \dontrun{
#' require(randomForest)
#' 
#' rfEstProc <- buildEstimationProcedure(train=randomForest)
#' rfEst <- rfEstProc(data=iris, .trainArgs=list(formula=Species ~ .),
#'            .predictArgs=list(type='prob'))
#'            
#' rfEst(iris[1:15, ])
#' rfEst(iris[1:15, ], list(type='response'))
#' 
#' form <- formula(Sepal.Length ~ Sepal.Width)
#' lmEstProc <- buildEstimationProcedure(lm)
#' lmEst <- lmEstProc(data=iris, .trainArgs = list(formula=form))
#' 
#' lmEst(newdata=NULL) # return fitted values
#' 
#' # reshaping GLM to predict outcomes in {-1,1}
#' glm_predict <- function(object, newdata) {
#'   2*round(predict(object, newdata, type='response')) - 1
#' }
#' 
#' Phi_glm <- buildEstimationProcedure(train=glm, predict=glm_predict)
#' }

buildEstimationProcedure <- function(train,
                                     predict=stats::predict,
                                     learningSet = "data",
                                     predictionSet = "newdata",
                                     modelName = "object") {
  f <- function(data, .trainArgs=NULL, .predictArgs=NULL) {

    .trainArgs <- c(.trainArgs, list(data=data))
    
    #update name of learningSet argument
    names(.trainArgs)[length(.trainArgs)] <- learningSet 
    
    model <- do.call(train,.trainArgs)
    
    predictArgs <- .predictArgs
    function(newdata, .predictArgs = predictArgs) {
      .predictArgs <- c(list(object=model, newdata=newdata), .predictArgs)
      names(.predictArgs)[1:2] <- c(modelName, predictionSet)
      
      do.call(predict, .predictArgs)
    }
  }
  class(f) <- c("estimationProcedure", class(f))
  f
}

