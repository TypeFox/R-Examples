jaatha_data_class <- R6::R6Class("jaatha_data", 
  lock_objects = TRUE, lock_class = TRUE,
  private = list(
    values = list(),
    options = list(),
    log_factorials = list()
  ),
  public = list(
    initialize = function(data, model) {
      assert_that(is_jaatha_model(model))
      private$values <- lapply(model$get_sum_stats(), function(stat) {
        private$options[[stat$get_name()]] <- stat$generate_data_opts(data)
        value <- stat$calculate(data, private$options[[stat$get_name()]])
        if (any(!is.finite(value))) {
          stop("Data has missing values for statistic ", 
               stat$get_name(), call. = FALSE)
        }
        value
      })
      private$log_factorials <- lapply(private$values, function(x) {
        log_factorials <- vapply(x, function(y) {
          if (y == 0) return(0)
          assert_that(is.count(y))
          sum(log(1:y))
        }, numeric(1))
        assert_that(all(is.finite(log_factorials)))
        log_factorials
      })
    },
    get_values = function(stat = NULL) {
      if (is.null(stat)) return(private$values)
      if (is.character(stat) || is.numeric(stat)) return(private$values[[stat]])
      private$values[[stat$get_name()]]
    },
    get_options = function(stat = NULL) {
      if (is.null(stat)) return(private$options)
      if (is.character(stat) || is.numeric(stat)) {
        return(private$options[[stat]])
      }
      private$options[[stat$get_name()]]
    },
    get_log_factorial = function(stat = NULL) {
      if (is.null(stat)) return(private$log_factorials)
      if (is.character(stat) || is.numeric(stat)) {
        return(private$log_factorials[[stat]])
      }
      private$log_factorials[[stat$get_name()]]
    },
    set_options = function(options) {
      warning("This function is intented for internal use only")
      private$options <- options
    }
  )
)


is_jaatha_data <- function(x) inherits(x, "jaatha_data")


#' Prepare the observed data for Jaatha
#' 
#' By default, this function assumes that the observed data is in a format 
#' identical to the format of the simulation results, before the summary
#' statistics are calculated.
#' 
#' @section Demographic Inference:
#' When used with \pkg{coala}, \code{coala::calc_sumstats_from_data()} should
#' create output that is compatible with this function.
#' 
#' @param data The data to be analysed with Jaatha. 
#'   It should be in a format identical to the 
#'   simulation results (see \code{\link{create_jaatha_model}}).
#' @param model The jaatha model, see \code{\link{create_jaatha_model}}.
#' @param ... Currently ignored.
#' @export
create_jaatha_data <- function(data, model, ...) UseMethod("create_jaatha_data")


#' @describeIn create_jaatha_data The data's format is identicial to the 
#'   simulated data.
#' @export
create_jaatha_data.default <- function(data, model, ...) {
  jaatha_data_class$new(data, model)
}


create_test_data <- function(model) {
  test_data <- model$test(quiet = TRUE)
  create_jaatha_data(test_data, model)
}
