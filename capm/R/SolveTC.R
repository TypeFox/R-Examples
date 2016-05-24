#' Modelling of reversible contraception for companion animals
#' @description System of ordinary differential equations to simulate the effect of reversible contraception in a population at equilibrium, where deaths are compensated by births and net immigration.
#' @param pars a named \code{\link{vector}} of length 5. The values are point estimates of the death rate (d), the fertility recovery rate (f), the sterilization rate (s), the proportion of infertile immigrants (z) and the proportion of the death rate compensated by immigration (r). Abreviations in parentheses indicate the names that must be given to the values.
#' @param init a named \code{\link{vector}} of length 2, with the total number of fertil (n) and infertil (g) animals.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param f.range optional sequence (between 0 and 1) with the fertility recovery rates to be simulated.
#' @param s.range optional \code{\link{vector}} of length 2, with a range of sterilization rates to be assessed. If given, the rates evaluated are those specified by the argument plus the point estimate given in \code{pars}.
#' @param z.range optional \code{\link{vector}} of length 2, with a range of the proportion of infertile immigrants. If given, the rates evaluated are those specified by the argument plus the point estimate given in \code{pars}.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @return \code{\link{list}}. The first element, \code{name}, is a string with the name of the function, the second element, \code{model}, is the model function. The third, fourth and fifth elements are vectors (\code{pars}, \code{init}, \code{time}, respectively) containing the \code{pars}, \code{init} and \code{time} arguments of the function. The sisxthth element \code{results} is a \code{\link{data.frame}} with up to as many rows as elements in time. The first fourth columns contain the time and the variables: \code{n}, \code{g} and \code{u}. When *.range arguments are given, additional columns contain the variables \code{f}, \code{s} and \code{z}.
#' @references \url{http://oswaldosantos.github.io/capm}
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' pars.solvetc <- c(d = 1 / 6, f = 0.5, s = 0.2, 
#'                    z = 0.2, r = 0.8)
#'
#' init.solvetc <- c(n = 950, g = 50)
#' 
#' # Solve for point estimates.
#' solvetc.pt <- SolveTC(pars = pars.solvetc, 
#'                       init = init.solvetc, 
#'                       time = 0:10, method = 'rk4')
#' 
#' # Solve for parameter ranges.
#' solvetc.rg <- SolveTC(pars = pars.solvetc, 
#'                       init = init.solvetc, 
#'                       time = 0:15,
#'                       f.range = seq(0, 1, 0.1), 
#'                       s.range = c(0.05, 0.4), 
#'                       z.range = c(0.05, 0.4),
#'                       method = 'rk4')
#'
SolveTC <- function(pars = NULL, init = NULL, time = NULL, f.range = NULL, s.range = NULL, z.range = NULL, ...) {
  if(!setequal(names(pars), c('d', 'f', 'r', 's', 'z', 'd'))) {
    stop('Values in pars must have the following names:\nd, f, s, z, d.')
  }
  if(!setequal(names(init), c('n', 'g'))) {
    stop('Values in init must have the following names:\nn, g.')
  }
  SolveTCfu <- function(pars = NULL, init = NULL, time = NULL) {
    init['u'] <- 0
    SolveTC.fu <- function(time, init, pars) {
      with(as.list(c(init, pars)), { 
        dn <- r * d * (n + g) + (1 - r) * d * (1 - z) * (n + g) + f * g -
          (d + s) * n
        dg <- (1 - r) * d * z * (n + g) + s * n - (d + f) * g
        du <- s * n    
        list(c(dn, dg, du))
      })
    }
    init <- c(init['n'], init['g'], init['u'])
    
    SolveTC.out <- ode(times = time, func = SolveTC.fu, 
                       y = init, parms = pars, ...)
    return(as.data.frame(SolveTC.out))
  }
  
  if (is.null(f.range) & is.null(s.range) & is.null(z.range)) {
    output <- SolveTCfu(pars = pars, init = init, time = time)
    SolveTC <- list(
      name = 'SolveTC',
      model = SolveTCfu,
      pars = pars,
      init = init,
      time = time,
      results = output)
    class(SolveTC) <- 'capmModels'
    return(SolveTC)
  } else {
    output <- NULL
    paras <- pars
    z.range <- c(z.range[1], pars['z'], z.range[2])
    s.range <- c(s.range[1], pars['s'], s.range[2])
    for (i3 in 1:length(z.range)) {
      for (i2 in 1:length(s.range)) {
        for (i1 in 1:length(f.range)) {
          paras[c('z', 's', 'f') ] <- 
            c(z.range[i3], s.range[i2], f.range[i1])
          output <- rbind(output,
                          SolveTCfu(pars = paras, 
                                    init = init, 
                                    time = time))
        }
      }
    }
    output <- data.frame(
      output,
      f = rep(f.range, each = length(time)),
      s = rep(s.range, 
                 each = length(time) * length(f.range)),
      z = rep(z.range, 
                   each = length(time) * length(f.range) * 
                     length(s.range)))
    SolveTC <- list(
      name = 'SolveTC',
      model = SolveTCfu,
      pars = pars,
      init = init,
      time = time,
      results = output)
    class(SolveTC) <- 'capmModels'
    return(SolveTC)
  }
}