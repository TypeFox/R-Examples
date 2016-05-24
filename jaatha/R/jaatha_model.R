#' @importFrom parallel mclapply
#' @importFrom assertthat is.error
jaatha_model_class <- R6::R6Class("jaatha_model", 
  lock_objects = FALSE, lock_class = TRUE,
  private = list(
    par_ranges = NA,
    sum_stats = list(),
    scaling_factor = 1,
    add_statistic = function(stat) {
      name <- stat$get_name()
      if (!is.null(private$sum_stats[[name]])) {
        stop("There already is a summary statistic with name '", name, 
             "' in the model")
      }
      private$sum_stats[[name]] <- stat
    },
    slow_sim = FALSE
  ),
  public = list(
    initialize = function(sim_func, par_ranges, sum_stats, 
                          scaling_factor, test) {
      assert_that(is.function(sim_func))
      if (length(formals(sim_func)) != 1) {
        stop("The simulation function must have exactly one argument.")
      }
      private$sim_func <- sim_func
      private$par_ranges <- par_ranges_class$new(par_ranges)
      assert_that(is.list(sum_stats))
      lapply(sum_stats, private$add_statistic)
      assert_that(is_single_numeric(scaling_factor))
      private$scaling_factor <- scaling_factor
      if (test) self$test()
      invisible(NULL)
    },
    simulate = function(pars, data, cores = 1) {
      "conducts a simulation for each parameter combination in pars"
      assert_that(is.matrix(pars))
      assert_that(ncol(pars) == private$par_ranges$get_par_number())
      assert_that(all(0 - 1e-5 <= pars & pars <= 1 + 1e-5))
      assert_that(is_jaatha_data(data))
      assert_that(is.count(cores))
      
      # Generate a seed for each simulation
      seeds <- sample_seed(length(pars) + 1)
      
      # Simulate
      sim_data <- mclapply(1:nrow(pars), function(i) {
        set.seed(seeds[i])
        sim_pars <- private$par_ranges$denormalize(pars[i, ])
        
        # Simulate and dump frames if an error occurs
        withCallingHandlers({
          sim_result <- private$sim_func(sim_pars)
        }, error = function(e) {
          error_dump <- tempfile("jaatha_frame_dump_", fileext = ".Rda")
          dump.frames("sim_error_dump")
          save("sim_error_dump", file = error_dump)
          stop(paste(e$message, "[Frame dump written to", error_dump, "]"), 
               call. = FALSE)
        })
        
        
        # Calculate Summary Statistics
        sum_stats <- lapply(private$sum_stats, function(sum_stat) {
          sum_stat$calculate(sim_result, data$get_options(sum_stat))
        })
        
        # Add the parameter values
        sum_stats$pars <- sim_pars
        sum_stats$pars_normal <- pars[i, ]
        
        sum_stats
      }, mc.preschedule = !private$slow_sim, mc.cores = cores)
      
      failed <- vapply(sim_data, is.error, logical(1))
      if (any(failed)) {
        lapply(which(failed), function(x) {
          warning("Simulation failed: ", sim_data[[x]])
        })
        stop("Simulations failed, check your model.")
      }
        
      set.seed(seeds[length(seeds)])
      sim_data
    },
    get_par_ranges = function() private$par_ranges,
    get_sum_stats = function() private$sum_stats,
    test = function(quiet = FALSE) {
      time <- system.time(
        sim_data <- private$sim_func(private$par_ranges$get_middle())
      )["elapsed"]
      
      if (!quiet) {
        if (time > 30) warning("Each simulation takes about ", round(time),
                               "s, Jaatha might run for a long time.")
        
        if (time < 1) message("A simulation takes less than a second")
        else message("A simulation takes about ", round(time), "s")
        
        if (time > 5) private$slow_sim <- TRUE
      }
      
      invisible(sim_data)
    },
    get_scaling_factor = function() private$scaling_factor,
    get_par_number = function() private$par_ranges$get_par_number(),
    get_sim_func = function() private$sim_func
  )
)


#' Specify a Model for a Jaatha Analysis 
#' 
#' This function can be used to create models for an analysis with Jaatha.
#' Models can be created using simulation function  
#' (see \code{\link{create_jaatha_model.function}}) or using a \pkg{coala} 
#' model (see \code{\link{create_jaatha_model.coalmodel}}).
#' 
#' @param x The primary argument. Can be a function used for simulations,
#'   or a coala model.
#' @param ... Additional parameters passed on to the dispatch function.
#' @param scaling_factor If your model is a down-scaled version of your data,
#'   you can indicated this using this value. The estimated expectation values
#'   are multiplied with this factor before the likelihood is calculated.
#' @param test A logical indicating whether a simulation is performed to test
#'   the model.
#' @export
create_jaatha_model <- function(x, ..., scaling_factor = 1, test = TRUE) {
  UseMethod("create_jaatha_model")
}


create_jaatha_model.default <- function(x, ..., 
                                        scaling_factor = 1, 
                                        test = TRUE) {
  stop("Can create a model from an object of class '", class(x), "'")
}


#' Specify a jaatha model using a simulation function
#' 
#' This is the usual way to specify a jaatha model. An detailed exampled on 
#' doing so is given in the `jaatha-intro` vignette.
#' 
#' @param x A simulation function. This function takes model parameters as 
#'   input, and returns the simulated data. The function must take exactly one 
#'   argument, which is a numeric vector of model parameters for which the 
#'   simulation should be conducted. The function should return the simulation
#'   results in an arbitrary format, that is then passed on to the summary
#'   statistics.
#' @param par_ranges A matrix stating the possible values for the model 
#'   parameters. The matrix must have one row for each parameter, and two
#'   columns which state the minimal and maximal possible value for the 
#'   parameter.
#' @param sum_stats A list of summary statistics created with 
#'   \code{\link{create_jaatha_stat}}. The simulation results will be passed
#'   to the statistics, which should convert them into a numeric vector.
#' @param ... Currently unused.
#' @inheritParams create_jaatha_model
#'
#' @export
#' @export create_jaatha_model.function
#' @examples 
#' create_jaatha_model(function(x) rpois(10, x),
#'                     par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2),
#'                     sum_stats = list(create_jaatha_stat("sum", sum)))
create_jaatha_model.function <- function(x, par_ranges, sum_stats, ...,
                                         scaling_factor = 1, 
                                         test = TRUE) {
  jaatha_model_class$new(x, par_ranges, sum_stats, 
                         scaling_factor = scaling_factor, test = test)
}


is_jaatha_model <- function(x) inherits(x, "jaatha_model")


create_test_model <- function() {
  create_jaatha_model(function(x) stats::rpois(10, x),
                      par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2),
                      sum_stats = list(stat_identity(), stat_sum()),
                      test = FALSE)
}
