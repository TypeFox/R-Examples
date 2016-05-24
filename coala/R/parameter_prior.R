prior_par_class <- R6Class("prior_par", inherit = named_par_class,
  private = list(prior = NA),
  public = list(
    initialize = function(name, prior) {
      super$initialize(name)
      if (!is.expression(prior)) stop("prior is not expression")
        private$prior <- prior
      },
      print = function() {
        cat(private$name, ": prior `",  as.character(private$prior), "`\n",
          sep = "")
      },
      sample = function() {
        eval(private$prior, envir = new.env())
      }
   )
)


#' @describeIn parameter Creates a named parameter with a prior
#'  distribution. Before each simulation, the expression for the prior
#'  is evaluated. The resulting value can be used in
#'  \code{\link{par_expr}} under the chosen name.
#'
#' @export
#' @param prior An expression. Evaluation this expression should give
#'   a sample from the prior distribution you want for the parameter.
#'   For example using \code{rnorm(1)} gives a standard normal prior.
par_prior <- function(name, prior) {
  prior_par_class$new(name, as.expression(substitute(prior)))
}


is.prior_par <- function(par) inherits(par, "prior_par")


sample_par_priors <- function(model) {
  par <- get_parameter(model)
  prior_par <- par[vapply(par, is.prior_par, logical(1))]
  if (length(prior_par) == 0) return(numeric(0))
  values <- vapply(prior_par, function(x) x$sample(), numeric(1))
  names(values) <- vapply(prior_par, function(x) x$get_name(), character(1))
  values
}
