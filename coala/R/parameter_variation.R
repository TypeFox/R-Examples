#' @importFrom R6 R6Class
variation_par_class <- R6Class("variation_par", inherit = parameter_class,
  private = list(
    base_par = list(),
    func = "variation",
    add_parameter = function(parameter) {
      if (is.numeric(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.character(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.par(parameter)) {
        private$base_par[[length(private$base_par) + 1]] <- parameter
        expr <- parameter$get_expression()
      } else {
        stop("Unexpected type of parameter")
      }
      expr
    }),
  public = list(
    initialize = function(parameter, variance) {
      expr_mean <- private$add_parameter(parameter)
      expr_var <- private$add_parameter(variance)
      private$expr <- parse(text = paste0(private$func, "(", expr_mean,
                                          ", ", expr_var, ")"))
    },
    get_base_par = function() private$base_par
  )
)

#' @importFrom stats rgamma
variation <- function(mean, variance) {
  rgamma(1, mean ^ 2 / variance, mean / variance)
}

is.par_variation <- function(object) inherits(object, "variation_par")


#' Variable Parameters
#'
#' This function can be used to let the values of a parameter vary between
#' the different loci. When used, the values for the enclosed parameter
#' will follow a gamma distribution with mean of the parameters original
#' value, and the variance specified as argument \code{variance}. This requires
#' that the original value is positive. When using this, the simulators
#' are called separately for each locus, which can dramatically increase the
#' time needed to simulate models with many loci.
#'
#' @param par A parameter whichs value will be made variable between the loci.
#' @param variance The variance of the gamma distribution, which the values used
#'   for simulation will follow.
#' @export
#' @seealso For parameters that are identical for all loci: \code{\link{parameter}}
#' @examples
#' model <- coal_model(5, 10) +
#'   feat_mutation(par_variation(par_const(5), 10)) +
#'   sumstat_nucleotide_div()
#' simulate(model)
par_variation <- function(par, variance) {
  variation_par_class$new(par, variance)
}
