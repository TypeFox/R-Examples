#' @rdname extractionFunctions 
#' @export
#' @title Extraction functions for boostr object attributes
#' @description
#' Access the various attributes of "\code{boostr}" objects through
#' these functions. See \code{\link{boostBackend}} for a description of every
#' \code{boostr} attribute.
#' @param boostrObj an object of class "\code{boostr}" -- most likely the output
#' of \code{\link{boost}} or \code{\link{boostBackend}}.
#' @return The attribute referenced to in the function's title. E.g., 
#' \code{extractEstimationProcedure} returns the stored estimation procedure.
#' \code{ensembleEstimators} returns the ensemble of estimators built during
#' \code{\link{boostBackend}}.
 ensembleEstimators <- function(boostrObj) { 
  attr(boostrObj, 'ensembleEstimators')
} 

#' @rdname extractionFunctions  
#' @export
 reweighterOutput <- function(boostrObj) {
  attr(boostrObj, 'reweighterOutput')
} 

#' @rdname extractionFunctions  
#' @export
extractPerformanceOnLearningSet <- function(boostrObj) {
  attr(boostrObj, 'performanceOnLearningSet')
}

#' @rdname extractionFunctions  
#' @export
extractCalcBoostrPerformance <- function(boostrObj) {
  attr(boostrObj, 'calcBoostrPerformance')
}

#' @rdname extractionFunctions  
#' @export
 estimatorPerformance <- function(boostrObj) { 
attr(boostrObj, 'estimatorPerformance')
} 

#' @rdname extractionFunctions  
#' @export
 oobVec <- function(boostrObj) { 
attr(boostrObj, 'oobVec')
} 

#' @rdname extractionFunctions  
#' @export
 extractReweighter <- function(boostrObj) { 
attr(boostrObj, 'reweighter')
} 

#' @rdname extractionFunctions  
#' @export
 reweighterArgs <- function(boostrObj) { 
attr(boostrObj, 'reweighterArgs')
} 

#' @rdname extractionFunctions  
#' @export
 extractAggregator <- function(boostrObj) { 
attr(boostrObj, 'aggregator')
} 

#' @rdname extractionFunctions  
#' @export
 aggregatorArgs <- function(boostrObj) { 
attr(boostrObj, 'aggregatorArgs')
} 

#' @rdname extractionFunctions  
#' @export
 extractEstimationProcedure <- function(boostrObj) { 
attr(boostrObj, 'estimationProcedure')
} 

#' @rdname extractionFunctions  
#' @export
 estimationProcedureArgs <- function(boostrObj) { 
attr(boostrObj, 'estimationProcedureArgs')
} 

#' @rdname extractionFunctions  
#' @export
 extractData <- function(boostrObj) { 
attr(boostrObj, 'data')
} 

#' @rdname extractionFunctions  
#' @export
 extractAnalyzePerformance <- function(boostrObj) { 
attr(boostrObj, 'analyzePerformance')
} 

#' @rdname extractionFunctions  
#' @export
 analyzePerformanceArgs <- function(boostrObj) { 
attr(boostrObj, 'analyzePerformanceArgs')
} 

#' @rdname extractionFunctions  
#' @export
extractSubsetFormula <- function(boostrObj) { 
  attr(boostrObj, 'subsetFormula')
} 

#' @rdname extractionFunctions  
#' @export
 extractFormatData <- function(boostrObj) { 
attr(boostrObj, 'formatData')
} 

#' @rdname extractionFunctions  
#' @export
 extractInitialWeights <- function(boostrObj) { 
attr(boostrObj, 'initialWeights')
} 