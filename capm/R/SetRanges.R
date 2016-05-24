#' Parameter ranges for global sensitivity analysis
#' @description Define the minimum and maximum values for parameters whose global sensitivities are to be assesses with \code{\link{CalculateGlobalSens}} or \code{\link{sensRange}} functions.
#' @param pars the same \code{pars} vector used in one of the following functions: \code{\link{SolveSI}} or \code{\link{SolveIASA}}.
#' @param range scale factor to define the minimum and maximum for each parameter. The default is 0.1, which set the minimum and maximum as 10 percent lesser and greater than the \code{pars} values.
#' @return \code{\link{data.frame}} with the complete set of parameter ranges.
#' @references Soetaert K and Petzoldt T (2010). Inverse modelling, sensitivity and monte carlo analysis in R using package FME. Journal of Statistical Software, 33(3), pp. 1-28.
#' 
#' Reichert P and Kfinsch HR (2001). Practical identifiability analysis of large environmental simulation models. Water Resources Research, 37(4), pp. 1015-1030.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @seealso \code{\link{sensRange}} and \code{\link{SolveSI}}.
#' @export
#' @examples
#' 
#' data(psu.ssu)
#' data(survey.data)
#' 
#' #####################
#' ## SolveIASA model ##
#' #####################
#' 
#' # Parameters and initial conditions.
#' pars.solve.iasa = c(
#'    b1 = 21871, b2 = 4374,
#'    df1 = 0.104, dm1 = 0.098, df2 = 0.125, dm2 = 0.118,
#'    sf1 = 0.069, sf2 = 0.05, sm1 = 0.028, sm2 = 0.05,
#'    k1 = 98050, k2 = 8055, h1 = 1, h2 = 0.5,
#'    a = 0.054, alpha = 0.1, v = 0.2, z = 0.1)
#'
#' # Set ranges 10 % greater and lesser than the 
#' # point estimates.
#' rg.solve.iasa <- SetRanges(pars.solve.iasa)
#' 
SetRanges <- function(pars = NULL, range = 0.1) {
  par.ranges <- data.frame(min = c(pars) * (1 - range),
                           max = c(pars) * (1 + range))
  rownames(par.ranges) = names(pars)
  return(par.ranges)
}
