#' Cross tabulations of expected and realized population sizes.
#' 
#' Compute contingency tables of expected (i.e., estimated) and realized (i.e.,
#' simulated) population sizes. The expected values are obtained with the
#' Horvitz-Thompson estimator.
#' 
#' The contingency tables are computed with \code{\link{tableWt}}.
#' 
#' @name spTable
#' @param inp an object of class \code{\linkS4class{simPopObj}} containing
#' household survey and simulated population data.
#' @param select character; vector defining the columns in slots 'pop' and
#' 'sample' of argument 'input' that should be used for tabulation.
#' @return A list of class \code{"spTable"} with the following components:
#' \item{expected}{the contingency table estimated from the survey data.}
#' \item{realized}{the contingency table computed from the simulated population
#' data.}
#' @note Sampling weights are automatically used from the input object 'inp'!
#' @author Andreas Alfons and Bernhard Meindl
#' @seealso \code{\link{spMosaic}}, \code{\link{tableWt}}
#' @keywords dplot
#' @export
#' @examples
#' 
#' set.seed(1234)  # for reproducibility
#' data(eusilcS)   # load sample data
#' samp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
#'   strata="db040", weight="db090")
#' eusilcP <- simStructure(data=samp, method="direct", basicHHvars=c("age", "rb090"))
#' res <- spTable(eusilcP, select = c("age", "rb090"))
#' class(res)
#' res
spTable <- function(inp, select) {
  if ( class(inp) != "simPopObj") {
    stop("wrong input! Argument 'inp' must be of class 'simPopObj'!\n")
  }
  dataS <- inp@sample@data
  dataP <- inp@pop@data
  weights <- inp@sample@weight
  
  # initializations
  # prepare weights
  if ( is.null(weights) ) {
    n <- nrow(dataS)
    weights <- rep.int(nrow(dataP)/n, n)
  } else {
    weights <- dataS[[weights]]
  }
  # prepare data seta
  if ( !all(select %in% colnames(dataS)) ) {
    stop("not all variables in argument 'select' available in the survey data available in argument 'inp'!\n")
  }
  if ( !all(select %in% colnames(dataP)) ) {
    stop("not all variables in argument 'select' available in the population data available in argument 'inp'!\n")
  }    
  dataS <- dataS[, select, with=FALSE]
  dataP <- dataP[, select, with=FALSE]

  # compute weighted table from sample (expected)
  tabS <- tableWt(dataS, weights)
  # compute table for population (realized)
  tabP <- table(dataP)
  # create object and return result
  res <- list(expected=tabS, realized=tabP)
  class(res) <- "spTable"
  res
}

# methods for class "spTable"
as.array.spTable <- function(x, ...) {
  values <- c(as.integer(x$expected), as.integer(x$realized))
  d <- c(dim(x$expected), 2)
  dnNew <- list(which=c("expected", "realized"))
  dn <- c(dimnames(x$expected), dnNew)
  array(values, dim=d, dimnames=dn)
}

as.table.spTable <- function(x, ...) {
  tab <- as.array(x)
  class(tab) <- "table"
  tab
}

plot.spTable <- function(x, ...) spMosaic(x, ...)

print.spTable <- function(x, ...) {
  # expected (from sample)
  cat("Expected:\n")
  print(x$expected, ...)
  # realized (population)
  cat("\nRealized:\n")
  print(x$realized, ...)
}
