#' Global sensitivity analysis
#' @description Wraper for \code{\link{sensRange}} function, which calculates population size sensitivities, to parameters used in one of the following functions: \code{\link{SolveIASA}}, \code{\link{SolveSI}} or \code{\link{SolveTC}}.
#' @param model.out an output from one of the previous function or a \code{\link{list}} with equivalent structure.
#' @param ranges output from the \code{\link{SetRanges}} function, applied to the \code{pars} argument used in the function previously specified in \code{model.out}.
#' @param sensv string with the name of the the output variables for which the sensitivity needs to be estimated.
#' @param all logical. If \code{\link{FALSE}}, sensitivity ranges are calculated for each parameter. If \code{TRUE}, sensitivity ranges are calculated for the combination of all aparameters.
#' @details When \code{all} is equal to TRUE, \code{dist} argument in \code{\link{sensRange}} is defined as "latin" and when equal to \code{\link{FALSE}}, as "grid". The \code{num} argument in \code{\link{sensRange}} is defined as 100.
#' @return A \code{data.frame} (extended by \code{summary.sensRange} when \code{all == TRUE}) containing the parameter set and the corresponding values of the sensitivity output variables.
#' @references Soetaert K and Petzoldt T (2010). Inverse modelling, sensitivity and monte carlo analysis in R using package FME. Journal of Statistical Software, 33(3), pp. 1-28.
#' 
#' Reichert P and Kfinsch HR (2001). Practical identifiability analysis of large environmental simulation models. Water Resources Research, 37(4), pp.1015-1030.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @seealso \code{\link{sensRange}}.
#' @export
#' @examples
#' #####################
#' ## SolveIASA model ##
#' #####################
#' 
#' ## Parameters and intial conditions.
#' pars.solve.iasa = c(
#'    b1 = 21871, b2 = 4374,
#'    df1 = 0.104, dm1 = 0.098, df2 = 0.125, dm2 = 0.118,
#'    sf1 = 0.069, sf2 = 0.05, sm1 = 0.028, sm2 = 0.05,
#'    k1 = 98050, k2 = 8055, h1 = 1, h2 = 0.5,
#'    a = 0.054, alpha = 0.1, v = 0.2, z = 0.1)
#'    
#' init.solve.iasa = c(
#'    f1 = 33425, fs1 = 10865,
#'    m1 = 38039, ms1 = 6808,
#'    f2 = 3343, fs2 = 109,
#'    m2 = 3804, ms2 = 68)
#'    
#' 
#' # Solve for point estimates.
#' solve.iasa.pt <- SolveIASA(pars = pars.solve.iasa, 
#'                           init = init.solve.iasa, 
#'                           time = 0:15, method = 'rk4')
#' 
#' ## Set ranges 10 % greater and lesser than the
#' ## point estimates.
#' rg.solve.iasa <- SetRanges(pars = pars.solve.iasa)
#' 
#' ## Calculate golobal sensitivity of combined parameters.
#' ## To calculate global sensitivity to each parameter, set
#' ## all as FALSE.
#' glob.all.solve.iasa <- CalculateGlobalSens(
#'   model.out = solve.iasa.pt,
#'   ranges = rg.solve.iasa, 
#'   sensv = 'n2', all = TRUE)
#' 
CalculateGlobalSens <- function(model.out = NULL, ranges = NULL, sensv = NULL, all = FALSE) {
  if (!setequal(rownames(ranges), names(model.out$pars))) {
    stop('All parameters in ranges must be\ncontained in "pars" argument  of model.out.')
  }
  if (length(intersect(sensv, names(model.out$init))) == 0) {
    stop('All variables in sensv must be\ncontained in "init" argument  of model.out.')
  }
  if (all) {
    sens <- sensRange(func = model.out$model, 
                      parms = model.out$pars, 
                      init = model.out$init,
                      time = model.out$time,
                      parRange = ranges,
                      dist = "latin", 
                      sensvar = sensv, num = 100)
    return(summary(sens))
  } else {
    sens <- NULL
    for (i in 1:length(model.out$pars)) {
      tmp <- sensRange(func = model.out$model, 
                       parms = model.out$pars, 
                       init = model.out$init,
                       time = model.out$time,
                       parRange = ranges[i, ],
                       dist = "grid", 
                       sensvar = sensv, num = 100)
      sens <- rbind(sens, summary(tmp))
    }
    param <- rep(names(model.out$pars), 
                 each = length(tmp[-1]))
    sens <- cbind(sens, param)
    return(sens)
  }
}