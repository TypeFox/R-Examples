#' @include model.R
#' @importFrom R6 R6Class
feature_class <- R6Class("feature", inherit = model_part,
  private = list(
    parameter = list(),
    population = NULL,
    time = NULL,
    rate = NA,
    set_population = function(population, expected_length = 1) {
      assert_that(identical(population, "all") || is.numeric(population))
      assert_that(identical(population, "all") ||
                    length(population) == expected_length)
      private$population <- population
    },
    add_parameter = function(parameter, required = TRUE, add_par = TRUE) {
      expr <- NA
      if (is.numeric(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.character(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.par(parameter)) {
        idx <- as.character(length(private$parameter) + 1)
        private$parameter[[idx]] <- parameter
        expr <- parameter$get_expression()
      } else if (is.null(parameter) || is.na(parameter)) {
        if (required) stop("Missing value for required parameter")
        return(NA)
      } else stop("Unexpected type of parameter")
      if (!add_par) return(expr)
      paste0("par(", expr, ")")
    }
  ),
  public = list(
    initialize = function(rate, population, time) {
      if (!missing(rate)) private$rate <- private$add_parameter(rate)
      if (!missing(time)) private$time <- private$add_parameter(time)
      if (!missing(population)) private$set_population(population)
    },
    get_parameters = function() private$parameter,
    get_population = function() private$population,
    get_time = function() private$time,
    reset_parameters = function() private$parameter <- list(),
    print = function() {
      cat("Feature of type", private$feature_table$type[1], "\n")
    },
    get_call = function() private$call,
    get_rate = function() private$rate,
    check = function(model) invisible(NULL)
  )
)

is.feature <- function(feature) inherits(feature, "feature")

ignore_par <- function(feature, model) ""
