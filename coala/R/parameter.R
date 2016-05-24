# Base class for all parameters.
# Contains an expression that can be assigned to some part of a feature.
#' @include model.R
#' @importFrom R6 R6Class
parameter_class <- R6Class("parameter", inherit = model_part,
  private = list(expr = NA),
  public = list(
    initialize = function(expr) {
      if (!is.expression(expr)) stop("No expression provided: ", expr,
                                     " is of type ", is(expr))
      private$expr <- expr
    },
    eval = function(envir = parent.frame()) {
      eval(private$expr, envir = envir)
    },
    get_expression = function() private$expr
  )
)


is.par <- function(par) {
  any("parameter" == class(par))
}



# Base class for Model parameters.
# Model parameters have a name, and a value is assigned to a variable of that
# name for each simulation.
#' @importFrom R6 R6Class
named_par_class <- R6Class("named_par", inherit = parameter_class,
  private = list(name = NA),
  public = list(
    initialize = function(name) {
      if (!(is.character(name) & length(name) == 1))
        stop("The parameter name must be a character")

      super$initialize(parse(text = name))
      private$name <- name
    },
    get_name = function() private$name,
    print = function() {
      cat(private$name, ": Named parameter")
    },
    check_value = function(value) TRUE,
    generate_value = function(pars) {
      if (is.null(names(pars)) || !any(private$name == names(pars))) {
        stop("No value for parameter ", private$name, " found")
      }
      value <- pars[[private$name]]
      self$check_value(value)
      value
    }
  )
)


is.named_par <- function(par) {
  any("named_par" == class(par))
}


#' Model Parameters
#'
#' These functions add parameters to a model. Parameters can either
#' be used in features, or added directly to a model using the plus operator.
#' The value of parameters can be specified in the simulation command
#' (for \code{par_named} and \code{par_range}), sampled from a prior
#' distribution (\code{par_prior}) or can be derived from other parameters
#' (\code{par_expr}).
#'
#' @name parameter
#' @seealso For parameters that variate between the loci in a model:
#'   \code{\link{par_variation}}, \code{\link{par_zero_inflation}}
#' @author Paul Staab
#' @examples
#' # A parameter (here for the mutation rate) that is always
#' # equal to '5':
#' model_base <- coal_model(20, 1) +
#'   sumstat_nucleotide_div()
#'
#' model <- model_base +
#'   feat_mutation(par_const(5))
#' simulate(model)
#'
#' # With using a prior:
#' model <- model_base +
#'   feat_mutation(par_prior("theta", rnorm(1, 5, .1)))
#' simulate(model)
#'
#' # Using a named parater:
#' model <- model_base +
#'   feat_mutation(par_named("theta"))
#' simulate(model, pars = c(theta = 5))
#'
#' # or similarly a ranged parameter:
#' model <- model_base +
#'   feat_mutation(par_range("theta", 1, 10))
#' simulate(model, pars = c(theta = 5))
#'
#' # Expressions can be used to derive parameters from
#' # other parameters:
#' model <- model_base +
#'   par_named("theta_half") +
#'   feat_mutation(par_expr(theta_half * 2))
#' simulate(model, pars = c(theta_half = 2.5))
#'
#' model <- model_base +
#'   par_named("theta_log") +
#'   feat_mutation(par_expr(exp(theta_log)))
#' simulate(model, pars = c(theta_log = log(5)))
NULL

#' @describeIn parameter Creates a parameter with value determined by evaluating an expression.
#' @param expr An R expression.
#'  This allows to define a parameter using an R expression.
#'  It can contain other named parameters (e.g. \code{2 * a} will create an
#'  parameter that is twice the value of an existing parameter \code{a}).
#'  Make sure that the expression always evaluates
#'  to a valid parameter value (a single numeric in almost all cases).
#' @export
par_expr <- function(expr) {
  parameter_class$new(as.expression(substitute(expr)))
}


#' @describeIn parameter Creates an parameter that is equal to a fixed value.
#'   Different to \code{par_expr}, the value is evaluated on parameter creation.
#' @export
#' @param constant An R expression.
#'   The constant value of the parameter.
#'   Different to \code{expr}, the expression is evaluated immediately and
#'   can not depend on other named parameters.
par_const <- function(constant) {
  parameter_class$new(as.expression(constant))
}


#' @describeIn parameter Creates an parameter whose value is specified via the
#'   \code{pars} argument in \code{\link{simulate.coalmodel}}.
#' @export
#' @param name Character. The name of the parameter. Must be unique in a model.
par_named <- function(name) {
  named_par_class$new(name)
}


range_par_class <- R6Class("range_par", inherit = named_par_class,
  private = list(range = NA),
  public = list(
    initialize = function(lower, upper, name) {
      stopifnot(all(is.numeric(c(lower, upper))))
      stopifnot(length(lower) == 1)
      stopifnot(length(upper) == 1)
      stopifnot(lower < upper)

      super$initialize(name)
      private$range <- c(lower, upper)
    },
    get_range = function() private$range,
    print = function() {
      cat(private$name, ": range between", private$range[1],
                           "and", private$range[2], "\n")
    },
    check_value = function(value) {
      if ((!is.numeric(value)) || length(value) != 1) {
        stop("The value for ", private$name, " must be a single numeric value")
      }
      if (value < private$range[1] - 1e-5 || value > private$range[2] + 1e-5) {
        stop("The value for ", private$name, " is not in the given range.")
      }
      TRUE
    }
  )
)


is.ranged_par <- function(par) inherits(par, "range_par")


#' @describeIn parameter Creates an parameter that can take a range of possible
#'  values.
#'  Similar to \code{\link{par_named}}, the value of the parameter
#'  used in a simulation is set via the \code{pars} argument.
#'  This is primarily intended for creating model parameters for
#'  \pkg{jaatha}.
#'
#' @export
#' @param lower A numeric. The lower boundary of the parameter's range.
#' @param upper A numeric. The upper boundary of the parameter's range.
par_range <- function(name, lower, upper) {
  range_par_class$new(lower, upper, name)
}


par_eval_func <- function(x) format(x, scientific = FALSE)
par_print_func <- function(x) substitute(x)

create_par_env <- function(model, parameters, ..., for_cmd = FALSE) {
  par_env <- new.env()

  if (!for_cmd) {
    par_env[["par"]] <- par_eval_func

    for (par in get_parameter(model)) {
      par_env[[par$get_name()]] <- par$generate_value(parameters)
    }
  } else {
    par_env[["par"]] <- par_print_func
  }

  additional_pars <- list(...)
  for (i in seq(along = additional_pars)) {
    par_env[[names(additional_pars)[i]]] <- additional_pars[[i]]
  }

  par_env[["zero_frac"]] <- create_zero_frac_func(additional_pars$locus_id,
                                                  additional_pars$locus_number)

  par_env
}


prepare_pars <- function(pars, model) {
  assert_that(is.numeric(pars))

  # Try to add names for non-prior parameters if missing
  if (is.null(names(pars))) {
    par_names <- get_par_names(model, without_priors = TRUE)
    if (length(pars) != length(par_names)) {
      stop("Unexpected number of parameters")
    }
    names(pars) <- par_names
  }

  # Sample from priors and return
  c(pars, sample_par_priors(model))
}


print_par <- function(par) paste0("`", substr(par, 5, nchar(par) - 1), "`")
