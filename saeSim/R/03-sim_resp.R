#' Response component
#' 
#' One of the components which can be added to a \code{sim_setup}.
#' 
#' @inheritParams sim_agg
#' @param respFun a function constructing the response variable
#' @inheritParams dplyr::mutate
#' 
#' @details Potentially you can define an \code{respFun} yourself. Take care that it only has one argument, named \code{dat}, and returns the a \code{data.frame}.
#' 
#' @seealso \code{\link{agg_all}}, \code{\link{sim_gen}}, \code{\link{sim_comp_pop}}, \code{\link{sim_sample}}, , \code{\link{sim_comp_sample}}
#' 
#' @rdname sim_resp
#' @export
#' @examples
#' base_id() %>% sim_gen_x() %>% sim_gen_e() %>% sim_resp_eq(y = 100 + 2 * x + e)
sim_resp <- function(simSetup, respFun) {
  sim_setup(simSetup, new("sim_fun", order = 3, call = match.call(), respFun))
}

#' @rdname sim_resp
#' @export
sim_resp_eq <- function(simSetup, ...) {
  
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- quote(mutate_wrapper)
  mc[[2]] <- NULL
  
  sim_setup(simSetup, new("sim_fun", order = 3, call = match.call(), eval(mc)))
}