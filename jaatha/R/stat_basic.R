stat_basic_class <- R6::R6Class("jaatha_stat_basic",
  lock_objects = FALSE, lock_class = TRUE,
  public = list(
    initialize = function(name, calc_func) {
      assert_that(is.character(name) && length(name) == 1)
      private$name <- name
      
      assert_that(is.function(calc_func))
      self$calculate <- calc_func
    },
    get_name = function() private$name,
    generate_data_opts = function(data) NULL
  ),
  private = list(
    name = ""
  )
)


#' Create a summary statistic for Jaatha
#' 
#' This function creates summary statistics for Jaatha models. A summary
#' statistic consists primarily of a function that calculates the statistic
#' from the simulation results. Jaatha primarily supports Poisson distributed
#' summary statistics, but can also transform summary statistics that follow
#' a different distribution in approximately Poisson distributed statistics.
#' 
#' @section Transformation of non Poisson distributed statistics:
#'   To transform a statistic into approximately Poisson distributed values,
#'   we first calculate the empirical quantiles of the real data for the 
#'   probabilities given in \code{breaks}. These are used as break points for
#'   divining the range of the statistic into disjunct intervals. We then count
#'   who many of the values for the simulated data fall into each intervals, and
#'   use this counts as summary statistic. The counts are multinomial 
#'   distributed, and should be close to the required Poisson distribution in 
#'   most cases.
#'     
#' @param name The name of the summary statistic
#' @param calc_func The function that summarizes the simulation data. Must take
#'   two arguments. The first is the simulated data, and the second are
#'   options that can be calculated from the real data. Ignoring the second 
#'   argument in the function body should be fine in most situations. The 
#'   function must return a numeric vector if \code{poisson = TRUE}, and can
#'   also return a numeric matrix if \code{poisson = FALSE}.
#' @param poisson If \code{TRUE}, it is assumed that the summary statistic
#'   values are (at least approximately) independent and Poisson distributed.
#'   If it is set to \code{FALSE}, the statistic is transformed into an approximately
#'   Poisson distributed array using a binning approach. See "Transformation
#'   of non Poisson distributed statistics" for details. If any summary
#'   statistic is only approximately Poisson distributed, Jaatha is a
#'   composite-likelihood method.
#' @param breaks The probabilities for the quantiles that are used for binning
#'   the data. See the section on non Poisson distributed summary statistics
#'   for details.
#' @return The summary statistic. Indented for being used with
#'   \code{\link{create_jaatha_model}}.
#' @export
create_jaatha_stat <- function(name, calc_func, poisson = TRUE, 
                               breaks = c(.1, .5, .9)) {
  
  if (poisson) return(stat_basic_class$new(name, calc_func))
  stat_cube_class$new(name, calc_func, breaks)
}


stat_identity <- function() create_jaatha_stat("id", function(x, y) x)
stat_sum <- function() create_jaatha_stat("sum", function(x, y) sum(x))


# @importFrom reshape2 melt
# Stat_PoiSmooth <- R6Class("Stat_PoiSmooth", inherit = stat_basic_class,
#   public = list(
#     get_model = function() private$model,
#     transform = function(sim_data) {
#       private$to_data_frame(sim_data$data)
#     },
#     initialize = function(data, name, model) {
#       super$initialize(data, name)
#       private$model = model
#     }
#   ),
#   private = list(
#     model = "",
#     to_data_frame = function(data) {
#       stopifnot(!is.null(data))
#       dim_names <- lapply(dim(data), function(x) 1:x)
#       names(dim_names) <- paste0("X", 1:length(dim(data)))
#       dimnames(data) <- dim_names
#       melt(data, value.name = "sum.stat")
#     }
#   )
# )
