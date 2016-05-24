stat_cube_class <- R6::R6Class("stat_cube", inherit = stat_basic_class,
  lock_objects = FALSE, lock_class = TRUE,
  private = list(
    break_values = numeric(),
    calc_break_values = function(values, props) {
      "calculate the actual values for break_values"
      break_values <- unique(quantile(values, props, na.rm = TRUE))
      break_values[break_values == 0] <- ifelse(length(break_values) == 1, 
                                    0.01, min(0.01, min(break_values[-1]) / 2))
      break_values
    },
    generate_cube = function(stat, break_values, cols) {
      assert_that(is.matrix(stat))
      assert_that(is.list(break_values))
      assert_that(ncol(stat) == length(break_values))
      assert_that(all(cols <= ncol(stat)))
      
      # Remove rows that contain NAs or NaNs
      stat <- stat[apply(stat, 1, function(x) all(is.finite(x))), , #nolint
                   drop = FALSE]
      
      # Classify the loci accordingly to their statistics
      locus_class <- matrix(1, nrow(stat), ncol(stat))
      for (i in 1:ncol(stat)) {
        for (brk in break_values[[i]]) {
          locus_class[, i] <- locus_class[, i] + (stat[, i] > brk)
        }
      }
      
      # Count the classes and return as vector
      dims <- vapply(break_values, length, numeric(1)) + 1
      factors <- cumprod(c(1, dims[-length(dims)]))
      classes_int <- apply(locus_class, 1, 
                           function(x) sum((x - 1) * factors) + 1) #nolint
      tabulate(classes_int, nbins = prod(dims))
    }
  ),
  public = list(
    initialize = function(name, calc_func, break_values) {
      assert_that(is.function(calc_func))
      private$calculate_matrix <- calc_func
      
      super$initialize(name, function(data, opts) {
        stat <- private$calculate_matrix(data)
        assert_that(is.numeric(stat))
        if (!is.matrix(stat)) stat <- matrix(stat, ncol = 1)
        private$generate_cube(stat, opts$break_values, 1:ncol(stat))
      })
      
      assert_that(is.numeric(break_values))
      assert_that(length(break_values) > 0)
      if (any(break_values < 0 | break_values > 1)) {
        stop("probs greater then one")
      }
      private$break_values <- break_values
    },
    generate_data_opts = function(data) {
      data_matrix <- private$calculate_matrix(data)
      assert_that(is.numeric(data_matrix))
      if (!is.matrix(data_matrix)) data_matrix <- matrix(data_matrix, ncol = 1)
      list(break_values = lapply(1:ncol(data_matrix), function(i) {
        private$calc_break_values(data_matrix[, i], private$break_values)
      }))
    }
  )
)
