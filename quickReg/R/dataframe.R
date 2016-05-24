#' Generic function for class 'reg'
#'
#' Generic function for class 'reg' to extracte concentrated results of regression models.
#' @param x A reg object
#' @param \dots additional arguments
#' @export
#' @return A data frame of univariate regression  result
#' @seealso \code{\link{dataframe.reg}},   \code{\link{reg}}
#' @examples
#' data(diabetes)
#' head(diabetes)
#'
#' reg_coxph<-reg(data = diabetes, x = c(1,4, 6),
#' y = 5, time = 2, factor = c(1, 3, 4), model = 'coxph')
#'
#' dataframe(reg_coxph,save=FALSE)

dataframe <- function(x, ...) UseMethod("dataframe") 
