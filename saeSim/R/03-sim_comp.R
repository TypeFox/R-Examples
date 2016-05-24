#' Calculation component
#' 
#' One of the components which can be added to a \code{sim_setup}. These functions can be used for adding new variables to the data.
#' 
#' @inheritParams sim_agg
#' @param fun a function, see details.
#' @param by names of variables as character; identifying groups for which fun is applied.
#' 
#' @details Potentially you can define a function for computation yourself. Take care that it only has one argument, named \code{dat}, and returns a \code{data.frame}. Use \code{\link{comp_var}} for simple data manipulation. Functions added with \code{sim_comp_pop} are applied before sampling; \code{sim_comp_sample} after sampling. Functions added with \code{sim_comp_agg} after aggregation.
#' 
#' @seealso \code{\link{comp_var}}, \code{\link{sim_gen}}, \code{\link{sim_agg}}, \code{\link{sim_sample}}, \code{\link{sim_comp_N}}, \code{\link{sim_comp_n}}, \code{\link{sim_comp_popMean}}, \code{\link{sim_comp_popVar}}
#' @export
#' @rdname sim_comp
#' @examples
#' # Standard behavior
#' sim_base() %>% sim_gen_x() %>% sim_comp_N()
#' 
#' # Custom data modifications
#' ## Add predicted values of a linear model
#' library(saeSim)
#'
#' comp_lm <- function(dat) {
#'   dat$linearPredictor <- predict(lm(y ~ x, data = dat))
#'   dat
#' }
#'
#' sim_base_lm() %>% sim_comp_pop(comp_lm)
#' 
#' # or if applied after sampling
#' sim_base_lm() %>% sim_sample() %>% sim_comp_pop(comp_lm)
sim_comp_pop <- function(simSetup, fun = comp_var(), by = "") {
  fun <- if(by == "") fun else apply_by(by, fun)
  sim_setup(simSetup, new("sim_fun", order = 4, call = match.call(), fun))
}

#' @export
#' @rdname sim_comp
sim_comp_sample <- function(simSetup, fun = comp_var(), by = "") {
  fun <- if(by == "") fun else apply_by(by, fun)
  sim_setup(simSetup, new("sim_fun", order = 6, call = match.call(), fun))
}

#' @export
#' @rdname sim_comp
sim_comp_agg <- function(simSetup, fun = comp_var(), by = "") {
  fun <- if(by == "") fun else apply_by(by, fun)
  sim_setup(simSetup, new("sim_fun", order = 8, call = match.call(), fun))
}
