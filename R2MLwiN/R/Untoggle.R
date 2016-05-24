#' Converts a categorical variable into several separate binary variables
#' 
#' This function converts a vector (factor) of categorical character strings
#' (integers) into several separate vectors of binary indicators to enable
#' back-compatibility with versions of \pkg{R2MLwiN} prior to 0.8-0.
#' 
#' @param categrv A vector (factor) of categorical character strings
#' (integers).
#' @param name A character string specifying a name of a prefix to be appended
#' the categories when generating dummy variables. If \code{NULL}, \code{v__} is
#' used as a prefix.
#' 
#' @return A matrix containing the generated dummy variables.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @examples
#' 
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: Assumes MLwiN path is C:/Program Files (x86)/MLwiN v2.30/
#' # ...so please change relevant line if different
#' # if using R2MLwiN via WINE, the path may look like 
#' # options(MLwiN_path='/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN v2.30/') 
#' 
#' # Example: tutorial
#' data(tutorial)
#' names(tutorial)
#' tutorial = cbind(tutorial, Untoggle(tutorial$school, 'school'))
#' names(tutorial)
#' }
#' 
#' @export
Untoggle <- function(categrv, name = NULL) {
  ## this function will untoggle categorical variable into a few separate binary variables
  vars <- unique(categrv)
  N <- length(vars)
  rvs <- sapply(1:N, function(x) as.integer(categrv == vars[x]))
  vals <- suppressWarnings(as.numeric(vars))
  cnames <- as.character(vars)
  if (!is.null(name)) {
    cnames[!is.na(vals)] <- paste(name, cnames[!is.na(vals)], sep = "_")
  } else {
    cnames[!is.na(vals)] <- paste("v_", cnames[!is.na(vals)], sep = "_")
  }
  colnames(rvs) <- cnames
  rvs
} 
