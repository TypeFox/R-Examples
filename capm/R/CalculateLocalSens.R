#' Local sensitivity analysis
#' @description Wraper for \code{\link{sensFun}} function, which estimates local effect of all model parameters on population size, applying the so-called sensitivity functions. The set of parameters used in any of the following functions can be assessed: \code{\link{SolveIASA}}, \code{\link{SolveSI}} or \code{\link{SolveTC}}.
#' @param model.out an output from one of the previous functions or a \code{\link{list}} with equivalent structure.
#' @param sensv string with the name of the output variables for which sensitivity needs to be estimated.
#' @details For further arguments of \code{\link{sensFun}}, defaults are used. See the help page of this function for details. Methods for class "sensFun" can be used.
#' @return a \code{\link{data.frame}} of class \code{\link{sensFun}} containing the sensitivity functions. There is one row for each sensitivity variable at each independent time. The first column \code{x}, contains the time value; the second column \code{var}, the name of the observed variable; and remaining columns have the sensitivity parameters.
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
#' ## Calculate local sensitivities to all parameters.
#' local.solve.iasa2 <- CalculateLocalSens(
#'   model.out = solve.iasa.pt, sensv = 'n2')
#' local.solve.iasa1 <- CalculateLocalSens(
#'   model.out = solve.iasa.pt, sensv = 'n1')
#' 
CalculateLocalSens <- function(model.out = NULL, sensv = 'n') {
  if (length(intersect(sensv, names(model.out$init))) == 0) {
    stop('All variables in sensv must be\ncontained in "init" argument  of model.out.')
  }
  sensFun(func = model.out$model, 
          parms = model.out$pars, 
          init = model.out$init,
          time = model.out$time,
          sensvar = sensv,
          varscale = 1)
}