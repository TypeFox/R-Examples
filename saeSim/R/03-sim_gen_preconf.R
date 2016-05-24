#' Preconfigured generation components
#' 
#' These are some preconfigured generation components and all wrappers around \code{\link{sim_gen}} and \code{\link{sim_gen_cont}}.
#' 
#' @inheritParams sim_gen_cont
#' @inheritParams sim_gen
#' @inheritParams sim_agg
#' @inheritParams gen_norm
#' 
#' @details \code{x}: fixed-effect component; \code{e}: model-error; \code{ec}: contaminated model error; \code{v}: random-effect (error constant for each domain); \code{vc} contaminated random-effect. Note that for contamination you are expected to add both, a non-contaminated component and a contaminated component.
#' @rdname sim_gen_preconf
#' @export
sim_gen_x <- function(simSetup, mean = 0, sd = 4, name = "x") {
  sim_gen(simSetup, generator = gen_norm(mean, sd, name))
}

#' @rdname sim_gen_preconf
#' @export
sim_gen_e <- function(simSetup, mean = 0, sd = 4, name = "e") {
  sim_gen(simSetup, generator = gen_norm(mean, sd, name))
}

#' @rdname sim_gen_preconf
#' @export
sim_gen_ec <- function(simSetup,
                       mean=0, sd=150, name = "e", 
                       nCont = 0.05, type = "unit", areaVar = "idD", fixed = TRUE) {
  sim_gen_cont(simSetup, generator = gen_norm(mean, sd, name), nCont = nCont, type = type, areaVar = areaVar, fixed = fixed)
}

#' @rdname sim_gen_preconf
#' @export
sim_gen_v <- function(simSetup, mean = 0, sd = 1, name = "v") {
  force(mean); force(sd); force(name)
  sim_gen(simSetup, gen_v_norm(mean, sd, name))
}

#' @rdname sim_gen_preconf
#' @export
sim_gen_vc <- function(simSetup,
                        mean=0, sd=40, name = "v", 
                        nCont = 0.05, type = "area", areaVar = "idD", fixed = TRUE) {
  force(mean); force(sd); force(name)
  sim_gen_cont(simSetup, generator = gen_v_norm(mean, sd, name), nCont = nCont, type = type, areaVar = areaVar, fixed = fixed)
}
