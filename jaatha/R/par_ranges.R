par_ranges_class <- R6::R6Class("par_ranges", 
  private = list(
    range = NA,
    log_range = NA,
    log_range_width = NA,
    offset = 0
  ),
  public = list(
    initialize = function(par_range) {
      assert_that(is.matrix(par_range))
      assert_that(ncol(par_range) == 2)
      assert_that(nrow(par_range) >= 1)
      assert_that(all(par_range[, 1] < par_range[, 2]))
      if (is.null(rownames(par_range))) {
        rownames(par_range) <- paste0("p", 1:nrow(par_range))
      }
      private$range <- par_range
      private$offset <- min(par_range) - 1
      private$log_range <- log(par_range - private$offset)
      private$log_range_width <- private$log_range[, 2] - 
                                   private$log_range[, 1]
    },
    normalize = function(value) {
      log_value <- log(value - private$offset)
      (log_value - private$log_range[, 1]) / private$log_range_width
    },
    denormalize = function(value) {
      exp(value * private$log_range_width + private$log_range[, 1]) + 
        private$offset
    },
    print = function() {
      print(private$range)
    },
    get_middle = function() {
      (private$range[, 2] - private$range[, 1]) / 2 + private$range[, 1]
    },
    get_par_number = function() nrow(private$range),
    get_par_names = function() rownames(private$range)
  )
)

is_par_ranges <- function(x) inherits(x, "par_ranges")
