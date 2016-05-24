#' @title Formula construction
#' @description Helper function for building multivariate scglr formula.
#' 
#' NOTE: Interactions involving factors are not allowed for now.
#' For interactions between two quantitative variables, use \code{I(x*y)} as usual.
#' @export
#' @param namesY a vector of character containing the names of the dependent variables.
#' @param namesX a vector of character containing the names of the covariates (X) involved
#'  in the components.
#' @param namesAX a vector of character containing the names of the additional covariates.
#' @return an object of class \code{Formula}.
multivariateFormula <- function(namesY,namesX,namesAX=NULL)
{
  form_lhs <- paste(namesY, collapse="+")
  form_rhs <- paste(namesX, collapse="+")
  if(!is.null(namesAX)) {
    form_rhs <- paste(form_rhs, paste(namesAX, collapse="+"), sep="|")
  }
  formula <- as.Formula(paste(form_lhs, "~", form_rhs, sep=""))
  environment(formula) <- .GlobalEnv
  return(formula)
}
