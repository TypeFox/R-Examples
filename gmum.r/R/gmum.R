# Silences R CMD check warnings
utils::globalVariables(c(".rs.getAnywhere.original", ".rs.getAnywhere"))
utils::globalVariables(c("X1", "X2", "label"))

#' @importFrom grDevices rainbow
#' @importFrom graphics hist pairs par plot plot.new points title
#' @importFrom stats na.omit prcomp predict rnorm runif update
#' @importFrom utils data read.csv
#' @importFrom Matrix Matrix
#' @importFrom methods as hasArg is new
#' @importFrom grid viewport pushViewport grid.newpage grid.layout 
NULL

# Lazy loading to allow for discovery of all files
evalqOnLoad( {
  # Autocompletion override
  autocompl <- function(x, pattern="") {
    targets <- c(asNamespace("Rcpp")$complete(x), x[['.staticFields']])
    grep(pattern, targets, value = TRUE)[! (substr(grep(pattern, targets, value = TRUE),1,1)==".")]
  }
  
  `.DollarNames.Rcpp_C++Object` <<- autocompl
  .DollarNames.Rcpp_SVMClient <<- autocompl
  .DollarNames.Rcpp_GNGServer <<- autocompl
  .DollarNames.Rcpp_CecModel <<- autocompl
  
  if(!exists(".rs.getAnywhere")) {
    .rs.getAnywhere <- NULL # Silences R CMD check warning
  }
  
  # Workaround RStudio bug
  if(exists(".rs.getAnywhere") && !exists(".rs.getAnywhere.original")) {
    .rs.getAnywhere.original <<- .rs.getAnywhere
    .rs.getAnywhere <<- function(a, envir=.GlobalEnv){ .rs.getAnywhere.original(a, .GlobalEnv) }
  }
})
