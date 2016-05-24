#' A Fuzzy logic evolutionnary algorithm.
#'
#' \tabular{ll}{
#' Package: \tab fugeR\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2012-07-11\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' This package allow you to use a genetic algorithm
#' in order to find a fuzzy system that
#' can be used as a prediction model.
#'
#' \code{\link{fugeR.run}} find a fuzzy system. \code{\link{fugeR.predict}}
#' make the prediction for the given input data.
#'
#' @name fugeR-package
#' @aliases fugeR-package
#' @docType package
#' @title A Fuzzy logic evolutionnary algorithm.
#' @author Alexandre Bujard \email{alexandre.bujard@@gmail.com}
#' @references
#' \url{http://library.epfl.ch/en/theses/?nr=2634}
#' @keywords package machine learning fuzzy logic system genetic algorithm model prediction expert decision data mining
#' @seealso \code{\link{fugeR.run}}
#' @examples
#'   \dontrun{
#'      fis <- fugeR.run (
#'                  In,
#'                  Out,
#'                  generation=100,
#'                  population=200,
#'                  elitism=40,
#'                  verbose=TRUE,
#'                  threshold=0.5,
#'                  sensiW=1.0,
#'                  speciW=1.0,
#'                  accuW=0.0,
#'                  rmseW=1.0,
#'                  maxRules=10,
#'                  maxVarPerRule=2,
#'                  labelsMf=2
#'              )
#'      
#'      prediction <- fugeR.predict( fis, In )
#'   }
NA
