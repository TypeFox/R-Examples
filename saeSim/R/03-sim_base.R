#' Base component
#' 
#' Use the `sim_base` functions to start a new \code{sim_setup}.
#' 
#' @param data a \code{data.frame}.
#' 
#' @export
#' @rdname sim_base
#' 
#' @examples
#' # Example for a linear model:
#' sim_base() %>% sim_gen_x() %>% sim_gen_e()
sim_base <- function(data = base_id(100, 100)) {
  new("sim_setup", base = data, simName = "")
}

#' Preconfigured set-ups
#' 
#' \code{sim_base_lm()} will start a linear model: One regressor, one error component. \code{sim_base_lmm()} will start a linear mixed model: One regressor, one error component and one random effect for the domain. \code{sim_base_lmc()} and \code{sim_base_lmmc()} add outlier contamination to the scenarios. Use these as a quick start, then you probably want to configure your own scenario.
#' 
#' @details Additional information on the generated variables:
#' \itemize{
#' \item{nDomains: }{100 domains}
#' \item{nUnits: }{100 in each domain}
#' \item{x: }{is normally distributed with mean of 0 and sd of 4}
#' \item{e: }{is normally distributed with mean of 0 and sd of 4}
#' \item{v: }{is normally distributed with mean of 0 and sd of 1, it is a constant within domains}
#' \item{e-cont: }{as e; probability of unit to be contaminated is 0.05; sd is then 150}
#' \item{v-cont: }{as v; probability of area to be contaminated is 0.05; sd is then 40}
#' \item{y }{= 100 + x + v + e}
#' }
#' @rdname sim_setup_preconfigured
#' @export
#' 
#' @examples
#' # The preconfigured set-ups:
#' sim_base_lm()
#' sim_base_lmm()
#' sim_base_lmc()
#' sim_base_lmmc()
sim_base_lm <- function() {
  sim_base() %>% 
    sim_gen_x(0, 4, name = "x") %>% 
    sim_gen_e(0, 4, name = "e") %>%
    sim_resp(function(dat) mutate_(dat, y = "100 + x + e"))
}

#' @rdname sim_setup_preconfigured
#' @export
sim_base_lmm <- function() {
  sim_base_lm() %>% sim_gen_v(0, 1, name = "v") %>% 
    sim_resp(function(dat) mutate_(dat, y = "y + v"))
}

#' @rdname sim_setup_preconfigured
#' @export
sim_base_lmc <- function() {
  sim_base_lm() %>% sim_gen_ec()
}

#' @rdname sim_setup_preconfigured
#' @export
sim_base_lmmc <- function() {
  sim_base_lmm() %>% 
    sim_gen_ec() %>%
    sim_gen_vc()
}
