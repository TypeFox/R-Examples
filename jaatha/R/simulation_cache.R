sim_cache_class <- R6::R6Class("sim_cache",
  private = list(
    max_size = 0,
    size = 0,
    sim_data = list(),
    deactivated = FALSE
  ),
  public = list(
    initialize = function(max_size = 10000) {
      assert_that(is.count(max_size) || max_size == 0)
      if (max_size == 0) private$deactivated <- TRUE
      private$max_size <- max_size
    },
    add = function(sim_data) {
      "Adds the simulations in sim_data to the cache"
      assert_that(is.list(sim_data))
      if (private$deactivated) {
        private$sim_data <- sim_data
        return(invisible(NULL))
      }
      
      is_valid <- vapply(sim_data, function(x) {
        if (!is.list(x)) return(FALSE)
        if (!any(names(x) == "pars_normal")) return(FALSE)
        TRUE
      }, logical(1))
      if (!all(is_valid)) stop("Invalid simulation data")
      
      private$size <- min(private$size + length(sim_data), private$max_size)
      private$sim_data <- c(sim_data, private$sim_data)[1:private$size]
      invisible(NULL)
    },
    get_sim_data = function(block) {
      "Returns all simulation that were conducted with the given block"
      if (private$deactivated) return(private$sim_data)
      
      in_block <- vapply(private$sim_data, function(x) {
        block$includes(x$pars_normal)
      }, logical(1))
      private$sim_data[in_block]
    },
    get_size = function() private$size,
    clear = function() {
      private$size <- 0
      private$sim_data <- list()
    }
  )
)

create_sim_cache <- sim_cache_class$new
