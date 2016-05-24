#'
#' @title Generic k-fold Cross Validation wrapper
#' 
#' @description A general abstraction of the k-fold 
#' cross validation procedure.
#' 
#' @usage
#' kFoldCV(proc, k, data, params,
#'        .rngSeed = 1234, .chunkSize = 1L, .doSEQ = FALSE)
#' 
#' @param proc the procedure to be k-fold cross validated. \code{proc} needs to
#'  accept \code{data} and \code{newdata} in its signature, and must return a
#'  numeric vector.
#' @param k the number of folds.
#' @param data a matrix or data.frame from which the folds will be created.
#' @param params a list or data.frame. If \code{params} is a list, every
#' combination of the entries in its cells will be used as parameters to be
#' cross validated. If \code{params} is a data.frame, each row of 
#' arguments will be cross-validated.
#' @param .rngSeed the seed set before randomly generating fold indices.
#' @param .chunkSize the number of parameter combinations to be processed
#' at once (see help for \code{iter}).
#' @param .doSEQ logical flag indicating whether cross validation should 
#' be run sequentially or with \code{\%dopar\%}. 
#' 
#' @details This function leverages \code{\link{foreach}} and 
#' \code{iter} to perform \code{k}-fold cross  validation in a distributed
#'  fashion (provided a parallel backend is registered).
#' 
#' Because the heart of this function is a pair of nested \code{foreach} loops
#' one should be careful of "over-parallelization". Meaning, if the routine
#' inside \code{proc} is already natively parallel, then by invoking this 
#' routine around \code{proc} you'll be distributing a distributed computation.
#' This may not yield the speed gains you would expect.
#' 
#' One work around to this -- assuming \code{proc} is parallelized using 
#' \code{foreach} is to call create a wrapper around \code{proc} that calls
#' \code{\link{registerDoSEQ}}. For example,
#'  
#' \code{proC <- function(...) {registerDoSEQ(); proc(...)}}
#' 
#' Alternatively, you could run \code{kFoldCV} sequentially by setting
#' \code{.doSEQ} to \code{TRUE}. 
#' 
#' For a procedure \code{proc <- function(data, newdata, arg1, ..., argN){...}}
#' , it may end up that cross-validating a single N-tuple of arguments
#' \code{c(arg1, ..., argN)} may be very quick. Hence, the time it takes 
#' to send off \code{proc}, the \code{data} and the appropriate combinations of
#' \code{params} may overwhelm the actual computation time. In this instance,
#' one should consider changing \code{.chunkSize} from 1 to \code{n}
#' (where \code{n} is any reasonable integer value that would justify the
#' passing of data to a distant node). 
#' 
#' @note
#' The current implementation of this assumes that entries in \code{params} are
#' numeric so that \code{as.matrix(expand.grid(params))} is a numeric matrix
#' with named columns. A work around to passing character parameters would be
#' to translate the character parameter to an integer, and write a wrapper
#' for \code{proc} that translates the interger back to the appropriate 
#' string. See the example below.
#' 
#' @return a vector whose length is equal to \code{nrow(params)}, if
#' \code{params} is a data.frame, or the number of combinations of elements of
#' \code{params} if it's a list. The i-th component corresponds to the k-fold
#' cross-validated value of \code{proc} evaluated with parameters from the i-th
#' combination of \code{params}.
#'
#' @export
#' 
#' @examples
#' # simple example with k-NN where we can build our own wrapper
#' library(class)
#' data(iris)
#' .iris <- iris[, 5:1] # put response as first column
#' 
#' # make a wrapper for class::knn
#' f <- function(data, newdata, k) {
#'   preds <- knn(train=data[,-1],
#'                test=newdata[, -1],
#'                cl=data[, 1],
#'                k=k)
#'   mean(preds==newdata[, 1])
#' }
#' 
#' params <- list(k=c(1,3,5,7))
#' 
#' accuracy <- kFoldCV(f, 10, .iris, params, .rngSeed=407)
#' 
#' data.frame(expand.grid(params), accuracy=accuracy)
#' 
#' # look at a more complicated example:
#' # cross validate an svm with different kernels and different models
#' require(e1071)
#' g <- function(data, newdata, kernel, cost, gamma, formula) {
#'   kern <- switch(kernel, "linear", "radial", stop("invalid kernel"))
#'   form <- switch(formula,
#'                  as.formula(Species ~ .),
#'                  as.formula(Species ~ Petal.Length + Petal.Width),
#'                  as.formula(Petal.Length ~ .),
#'                  stop('invalid formula'))
#'   
#'    svmWrapper <- function(data, newdata, kernel, cost, gamma, form) {
#'                    svmObj <- svm(formula=form, data=data, kernel=kernel,
#'                                  cost=cost, gamma=gamma)
#'                    predict(svmObj, newdata)
#'                  }
#'   preds <- svmWrapper(data, newdata, kernel=kern, cost=cost,
#'                       gamma=gamma, form=form)
#'   
#'   if (formula != 3) {
#'     mean(preds == newdata[["Species"]])
#'   } else {
#'     mean((preds - newdata[["Petal.Length"]])^2)
#'   }
#' }
#' 
#' params <- list(kernel=1:2, cost=c(10,50), gamma=0.01, formula=1)
#' accuracy <- kFoldCV(g, 10, iris, params)
#' data.frame(expand.grid(params), metric=accuracy)
#' 
kFoldCV <- function(proc, k, data, params, .rngSeed = 1234,
                    .chunkSize = 1L, .doSEQ = FALSE) {
  # randomly partition data --
  # this will throw a warning if nrow(data) %% k != 0
  set.seed(.rngSeed)
  permutedIndices <- sample.int(nrow(data))
  folds <- matrix(permutedIndices, ncol=k)
  
  
  # build parameter matrix according to whether user passed in list of
  # parameters to search over, or prebuilt data.frame
  parameterMat <- switch(class(params),
                         list = as.matrix(expand.grid(params)),
                         data.frame = as.matrix(params),
                         stop("'params' must either be a list or data.frame."))
  
  `%op%` <- `%do%`
  
  if ( !.doSEQ && getDoParRegistered() ) `%op%` <- `%dopar%`
  
  pars <- NULL # instantiate local variable
  testFoldNum <- NULL # instantiate local variable
  
  foreach(pars=iter(parameterMat, by="row", chunksize=.chunkSize),
          .combine=cbind, .final=colMeans) %:%
    foreach(testFoldNum=seq.int(k), .combine=rbind) %op% {
              
              testFold <- folds[, testFoldNum]
              train <- data[-testFold, ]
              test <- data[testFold, ]

              tmp <- foreach(par=iter(pars, by="row")) %do%
                      do.call(proc,
                              c(list(data=train), list(newdata=test), par))
              
              do.call(cbind, tmp)
            }
}