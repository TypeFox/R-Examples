# helper function to search through .procArgs list for any
# entry named "formula" that boost will try and subset the data
# according.
findFormulaIn <- function(.procArgs) {
  flatList <- c(.procArgs, recursive=TRUE)
  listNames <- names(flatList)
  formulaIndex <- grep(pattern=".*[.]{0,1}formula$", x=listNames, perl=TRUE)
  if (length(formulaIndex) == 1) {
    as.formula(flatList[[formulaIndex]])
  } else if (length(formulaIndex) == 0) {
    NULL
  } else {
    stop("multiple formulae in .procArgs...")
  }
}

#' @title Boost an estimation procedure with a reweighter and aggregator.
#' @description
#' Perform the Boost algorithm on \code{proc} with \code{reweighter} and
#' \code{aggregator} and monitor estimator performance with 
#' \code{analyzePerformance}.
#' 
#' 
#' @param B the number of iterations to run.
#' @param reweighter a boostr compatible reweighter function.
#' @param aggregator a boostr compatible aggregator function.
#' @param proc a boostr compatible estimation procedure.
#' @param data the learning set to pass to \code{proc}. \code{data} is assumed
#' to hold the response variable in its first column.
#' @param initialWeights a vector of weights used for the first iteration of
#' the ensemble building phase of Boost.
#' @param .procArgs a named list of arguments to pass to \code{proc} in
#' addition to \code{data}.
#' @param .reweighterArgs a named list of arguments to pass to 
#' \code{reweighter} in addition to \code{proc}, \code{data} and 
#' \code{weights}. These are generally initialization values for other 
#' parameters that govern the behaviour of \code{reweighter}.
#' @param .aggregatorArgs a named list of arguments to pass to 
#' \code{aggregator} in addition to the output from \code{reweighter}. 
#' @param .storeData a boolean indicating whether the data should be stored in
#' the returned \code{boostr} object under the attribute "\code{data}". 
#' @param .calcBoostrPerformance a boolean indicating whether 
#' \code{analyzePerformance} should be used to monitor the performance of the
#' returned \code{boostr} object on the learning set. A value of 
#' \code{seq.int(nrow(data))} will be passed to \code{analyzePerformance} as
#' the \code{oobObs} argument.
#' @param .subsetFormula a \code{formula} object indicating how \code{data} is
#' to be subsetted. A formula of like "Type ~ ." will rearrange the columns of
#' \code{data} such that \code{data[,1] == data$Type}. By default, this value
#' is taken to be the value of the \code{formula} entry in \code{.procArgs}.
#' If multiple entries have the substring "formula" in their names, the search
#' will throw an error and you're advised to manually set \code{.subsetFormula}.
#' @param .formatData a boolean indicating whether the data needs to be
#' reformatted via \code{.subsetFormula} such that the response variable is in
#' the first column and the remaining columns are all predictor variables. This
#' is defaulted to \code{!is.null(.subsetFormula)}. 
#' @param analyzePerformance a boostr compatible performance analyzer.
#' @param .analyzePerformanceArgs a named list arguments to pass to 
#' \code{analyzePerformance} in addition to \code{prediction}, \code{response},
#' and \code{oobPbs}.
#' 
#' @details
#' For the details behind this algorithm, check out the paper at 
#' \url{http://pollackphoto.net/misc/masters_thesis.pdf}
#' 
#' @references
#' Steven Pollack. (2014). Boost: a practical generalization of AdaBoost
#' (Master's Thesis). \url{http://pollackphoto.net/misc/masters_thesis.pdf}
#' 
#' @note
#' \code{\link{wrapReweighter}}, \code{\link{wrapAggregator}}, 
#' \code{\link{wrapPerformanceAnalyzer}}, \code{\link{wrapProcedure}}, and
#' \code{\link{buildEstimationProcedure}} are all Wrapper Generators
#' designed to allow user implemented functions inside the \code{boostBackend}.
#' These functions are intelligently called from inside \code{\link{boost}}.
#' Thus, to minimize any sources of frustration, the recommended use of 
#' \code{boostBackend} is through \code{\link{boost}}.
#' 
#' @return a "\code{boostr}" object. The returned closure is the output of
#' \code{aggregator} on the collection of estimators built during the iterative
#' phase of Boost. This is intended to be a new estimator, and hence accepts
#' the argument \code{newdata}. However, the estimator also has attributes
#' \item{ensembleEstimators}{An ordered list whose components are the trained
#'  estimators.}
#' \item{reweighterOutput}{An ordered list whose components are the output of
#' \code{reweighter} at each iteration.}
#' \item{performanceOnLearningSet}{The performance of the returned boostr object
#' on the learning set, as measure by \code{analyzePerformance}. This is only
#' calculated if \code{.calcBoostrPerformance=TRUE}}
#' \item{estimatorPerformance}{An ordered list whose components are the output
#' of \code{analyzePerformance} at each iteration.}
#' \item{oobVec}{A row-major matrix whose \eqn{ij}{ij}-th entry indicates if
#' observation \eqn{j}{j} was used to train estimator \eqn{i}{i}.}
#' \item{reweighter}{The reweighter function used.}
#' \item{reweighterArgs}{Any additional arguments passed to \code{boostBackend}
#' for \code{reweighter}.}
#' \item{aggregator}{The aggregator function used.}
#' \item{aggregatorArgs}{Any additional arguments passed to \code{boostBackend}
#' for \code{aggregator}.}
#' \item{estimationProcedure}{The estimation procedure used.}
#' \item{estimationProcedureArgs}{Any additional arguments passed to 
#' \code{boostBackend} for \code{proc}.}
#' \item{data}{The learning set. Only stored if \code{.storeData = TRUE}.}
#' \item{analyzePerformance}{The performance analyzer used.}
#' \item{analyzePerformanceArgs}{Any additional arguments passed to 
#' \code{boostBackend} for \code{analyzePerformance}.}
#' \item{subsetFormula}{The value of \code{.subsetFormula}.}
#' \item{formatData}{The value of \code{.formatData}.}
#' \item{storeData}{The value of \code{.storeData}.}
#' \item{calcBoostrPerformance}{The value of \code{.calcBoostrPerformance}}
#' \item{initialWeights}{The initial weights used.}
#'
#'  The attributes can be accessed through the appropropriate 
#'  \code{\link[=ensembleEstimators]{extraction function}}.
#' @examples
#' \dontrun{
#' df <- within(iris, {
#'               Setosa <- factor(2*as.numeric(Species == "setosa") - 1)
#'               Species <- NULL
#'              })
#' 
#' form <- formula(Setosa ~ . )
#' df <- model.frame(formula=form, data=df)
#' 
#' # demonstrate arc-fs algorithm using boostr convenience functions
#' 
#' glmArgs <- list(.trainArgs=list(formula=form, family="binomial"))
#' 
#' # format prediction to yield response in {-1,1} instead of {0,1}
#' glm_predict <- function(object, newdata) {
#'   2*round(predict(object, newdata, type='response')) - 1
#'   }
#'   
#' Phi_glm <- buildEstimationProcedure(train=glm, predict=glm_predict)
#' 
#' phi <- boostBackend(B=3, data=df,
#'                      reweighter=adaboostReweighter,
#'                      aggregator=adaboostAggregator,
#'                      proc=Phi_glm,
#'                      .procArgs=glmArgs)
#'}       
boostBackend <- 
  function(B, reweighter, aggregator, proc,
           data, initialWeights, .procArgs,
           analyzePerformance = defaultOOBPerformanceAnalysis,
           .reweighterArgs = NULL,
           .aggregatorArgs = NULL,
           .analyzePerformanceArgs = NULL,
           .subsetFormula = findFormulaIn(.procArgs),
           .formatData = !is.null(.subsetFormula),
           .storeData = FALSE,
           .calcBoostrPerformance=TRUE) {
  
  
  # format data and set default OOB performance analysis if not already set.
  if (.formatData) {
    data <- model.frame(formula=.subsetFormula, data=data)
  }
  
  if (is.null(analyzePerformance)) {
    warning("analyzePerformance is NULL")
    analyzePerformance <- function(...) {return(NULL)}
  }
  
  # initialize ensembleStats
  # this will hold all output from ensemble building phase of Boost
  ensembleStats <- NULL 
  
  reducer <- function(old,new) {
    # capture newest output in prevItersOutput
    prevItersOutput <<- new
    
    # fold old into new, and then update ensembleStats
    out <- lapply(names(new), function(name) {
      if (inherits(new[[name]], 'function') || inherits(new[[name]], 'list')) {
        c(old[[name]], list(new[[name]]))
      } else {
        rbind(old[[name]], new[[name]])
      }
    })
    
    names(out) <- names(new)
    ensembleStats <<- out
    ensembleStats
  }

  # initialize "prevItersOutput" to be equal weights and other
  # reweighter input. The point is that a legit reweighter should output the
  # data its next iteration needs.
  prevItersOutput <- c(list(weights = initialWeights), .reweighterArgs)
  
  n <- nrow(data)
  
  ensembleStats <- foreach(i=seq.int(B), .combine=reducer,
                           .init=ensembleStats) %do% {
    # build the estimator
    #analyze performance
    
    weights <- prevItersOutput$weights
    
    ibObs <- sample(n, replace=TRUE, prob=weights) #ib == in bag
    oobObs <- seq.int(n)[!(seq.int(n) %in% unique(ibObs))] #oob == out of bag
    
    oobVec <- rep.int(0, n)
    oobVec[oobObs] <- 1 # oobVec is 1 at entry i if entry i is out of bag
    
    estimator <- do.call(proc, c(.procArgs, list(data=data[ibObs, ])))
    
    ### do estimator performance analysis:
    preds <- estimator(data[,-1])
    truth <- data[, 1]
    
    performance <- do.call(analyzePerformance,
                           c(list(prediction=preds,
                                  response=truth,
                                  oobObs=oobObs),
                             .analyzePerformanceArgs)
                           )
    
    reweighterArgs <- c(list(prediction=preds, response=truth),
                        prevItersOutput)
    
    reweighterOutput <- do.call(reweighter, reweighterArgs)
                             
    c(estimators=estimator, list(performance=performance), reweighterOutput,
      list(oobVec=oobVec))
  }
  

  out <- do.call(aggregator, c(ensembleStats, .aggregatorArgs))
  class(out) <- "boostr"
  
  # calculate boosted estimators performance on learning set
  if (.calcBoostrPerformance) {
    attr(out, 'performanceOnLearningSet') <- 
      do.call(analyzePerformance, 
              c(list(prediction=out(data[,-1]),
                     response=data[,1],
                     oobObs=seq.int(n)),
                .analyzePerformanceArgs)
      )
  }

  attr(out, 'ensembleEstimators') <- ensembleStats$estimators
  
  attr(out, 'reweighterOutput') <- ensembleStats[!(names(ensembleStats) %in%
                                      c("estimators", "performance", "oobVec"))]
  
  attr(out, 'estimatorPerformance') <- ensembleStats$performance
  attr(out, 'oobVec') <- ensembleStats$oobVec
  
  attr(out, 'reweighter') <- reweighter
  attr(out, 'reweighterArgs') <- .reweighterArgs
  
  attr(out, 'aggregator') <- aggregator
  attr(out, 'aggregatorArgs') <- .aggregatorArgs
  
  attr(out, 'estimationProcedure') <- proc
  attr(out, 'estimationProcedureArgs') <- .procArgs
  
  attr(out, 'data') <- NULL
  if (.storeData) attr(out, 'data') <- data
  
  attr(out, 'analyzePerformance') <- analyzePerformance
  attr(out, 'analyzePerformanceArgs') <- .analyzePerformanceArgs
  
  attr(out, 'subsetFormula') <- .subsetFormula
  attr(out, 'formatData') <- .formatData
  attr(out, 'storeData') <- .storeData
  attr(out, 'calcBoostrPerformance') <- .calcBoostrPerformance
  
  attr(out, 'initialWeights') <- initialWeights
  
  return(out)
}
