#' Modelling of immigration, abandonment, sterilization and adoption of companion animals
#' @description System of ordinary differential equations to simulate the effect of immigration of owned dogs, abandonment, sterilization of owned and stray dogs and adoption, on population dynamics.
#' @param pars a named \code{\link{vector}} of length 21, with point estimates of model parameters (see details).
#' @param init a named \code{\link{vector}} of length 8, with point estimates of model parameters (see details).
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param s.range optional sequence (between 0 and 1) of the sterilization rates to be simulated.
#' @param a.range optional \code{\link{vector}} of length 2, with range (ie, confidence interval) of abandonment rates to be assessed. If given, the rates evaluated are those specified by the argument plus the point estimate given in \code{pars}.
#' @param alpha.range optional \code{\link{vector}} of length 2, with range (ie, confidence interval) of adoption rates to be assessed. If given, the rates evaluated are those specified by the argument plus the point estimate given in \code{pars}.
#' @param v.range optional \code{\link{vector}} of length 2, with range of values of immigration rates to be assessed. This must be expressed as a percentage of owned animals carrying capacity.
#' @param s.fm logical. If \code{TRUE}, s.range is used for females and males and if \code{FALSE}, it is only used for females (for males, the point estimate given in \code{pars} is used.)
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details The \code{pars} argument must contain named values, using the following conventions: \code{1}: owned animals; \code{2}: stray animals; \code{f}: females; \code{m}: males. Then:
#' 
#'  
#' \code{b1} and \code{b2}: number of births.
#' 
#' \code{df1}, \code{dm1}, \code{df2} and \code{dm2}: death rate.
#' 
#' \code{sf1}, \code{sm1}, \code{sf2} and \code{sm2}: sterilization rate.
#' 
#' \code{k1} and \code{k2}: carrying capacity.
#' 
#' \code{h1} and \code{h2}: mean harem size.
#' 
#' \code{a}: abandonment rate.
#' 
#' \code{alpha}: adoption rate.
#' 
#' \code{v}: immigration rate.
#' 
#' \code{z}: proportion of sterilized immigrants.
#' 
#' The \code{init} argument must contain named values for the inital number of animals, using the following conventions: \code{1}: owned animals; \code{2}: stray animals; \code{f}: females; \code{m}: males; and \code{s}: sterilized. Then, the names are:
#' 
#' \code{f1}, \code{fs1}, \code{m1}, \code{ms1}, \code{f2}, \code{fs2}, \code{m2} and \code{ms2}.
#' 
#' If any range is specified (e.g \code{s.range}), the remaining ranges must be specified too (\code{a.range}, \code{alpha.range} and \code{v.range}).
#' The function is a wrapper around the defaults of \link[deSolve]{ode} function, whose help page must be consulted for details. An exception is the method argument, which here has "rk4" as a default.
#' @return \code{\link{list}}. The first element, \code{name}, is a string with the name of the function, the second element, \code{model}, is the model function. The third, fourth and fifth elements are vectors (\code{pars}, \code{init}, \code{time}, respectively) containing the \code{pars}, \code{init} and \code{time} arguments of the function. The sisxthth element \code{results} is a \code{\link{data.frame}} with up to as many rows as elements in time. The first column contain the time and subsequent columns contain the size of specific subpopulations, named according to conventions above. The \code{group} column differentiate between owned and strays. When *.range arguments are given, the last fourth columsn specify their instances.
#' @note Logistic growth models are not intended for scenarios in which
#' population size is greater than carrying capacity and growth rate is negative.
#' @references \url{http://oswaldosantos.github.io/capm}
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
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
#'                           time = 0:8, method = 'rk4')
#' 
#' # Solve for parameter ranges.
#' solve.iasa.rg <- SolveIASA(pars = pars.solve.iasa, 
#'                           init = init.solve.iasa, 
#'                           time = 0:8,
#'                           s.range = seq(0, .4, l = 15), 
#'                           a.range = c(0, .2), 
#'                           alpha.range = c(0, .2),
#'                           v.range = c(0, .1),
#'                           method = 'rk4')
#'
SolveIASA <- function(pars = NULL, init = NULL, time = NULL, s.range = NULL, a.range = NULL, alpha.range = NULL, v.range = NULL, s.fm = TRUE, ...) {
  if(!setequal(names(pars), c('b1', 'b2', 'df1', 'dm1', 
                              'df2', 'dm2', 'sf1', 'sf2', 
                              'sm1', 'sm2', 'k1', 'k2', 'h1',
                              'h2', 'a', 'alpha', 'v', 'z'))) {
    stop('Values in pars must have the following names:\nb1, b2, df1, dm1, df2, dm2, sf1, sf2, sm1, sm2, k1, k2, h1, h2, a, alpha, v, z')
  }
  if(!setequal(names(init), c('f1', 'fs1', 'm1', 'ms1',
                              'f2', 'fs2', 'm2', 'ms2'))) {
    stop('Values in init must have the following names:\nf1, fs1, m1, ms1, f2, fs2, m2, ms2')
  }
  init['n1'] <- sum(init[c('f1', 'm1')])
  init['ns1'] <- sum(init[c('fs1', 'ms1')])
  init['n2'] <- sum(init[c('f2', 'm2')])
  init['ns2'] <- sum(init[c('fs2', 'ms2')])
  init['N1'] <- sum(init[c('n1', 'ns1')])
  init['N2'] <- sum(init[c('n2', 'ns2')])
  init['N'] <- sum(init[c('N1', 'N2')])
  
  SolveIASAfu <- function(pars, init, time) {
    SolveIASA.fu <- function(time, init, pars) {
      with(as.list(c(init, pars)), {
        
        if (f1 + fs1 + m1 + ms1 <= k1) {
          omega1 <- f1 + fs1 + m1 + ms1
        } else {
          omega1 <- k1
        }
        
        if (f2 + fs2 + m2 + ms2 <= k2) {
          omega2 <- f2 + fs2 + m2 + ms2
        } else {
          omega2 <- k2
        }
        
        x1 <- (b1 * (h1 * m1 + f1)) / (2 * h1 * f1 * m1)
        wf1 <- (x1 * m1) / (m1 + f1 * h1 ^ (-1))
        w.f1 <- wf1 - (wf1 - df1) * (omega1 / k1)
        c.f1 <- df1
        q <- k1 * v * (1 - z) / 2
        qs <- k1 * v * z / 2
        
        d.f1 <- (w.f1 - c.f1 - sf1 - a) * f1 +
          (alpha * f2 + q) * (1 - (omega1 / k1))
        
        d.fs1 <- - (c.f1 + a) * fs1 + sf1 * f1 + 
          (alpha * fs2 + qs) * (1 - (omega1 / k1))
        
        wm1 <- (x1 * f1) / (m1 + f1 * h1 ^ (-1))
        w.m1 <- wm1 - (wm1 - dm1) * (omega1 / k1)
        c.m1 <- dm1
        
        d.m1 <- (w.m1 - c.m1 - sm1 - a) * m1 +
          (alpha * m2 + q) * (1 - (omega1 / k1))
        
        d.ms1 <- - (c.m1 + a) * ms1 + sm1 * m1 + 
          (alpha * ms2 + qs) * (1 - (omega1 / k1))
        
        x2 <- (b2 * (h2 * m2 + f2)) / (2 * h2 * f2 * m2)
        w.f2 <- (m2 * x2) /  (m2 + f2 * h2 ^ (-1))
        c.f2 <- df2 + (w.f2 - df2) * (omega2 / k2)
        
        d.f2 <- (w.f2 - c.f2 - sf2 - alpha) * f2 +
          a * f1 * (1 - (omega2 / k2))
        
        d.fs2 <- - (c.f2 + alpha) * fs2 + sf2 * f2 +
          a * fs1 * (1 - (omega2 / k2))
        
        w.m2 <- (f2 * x2) / (m2 + f2 * h2 ^ (-1))
        c.m2 <- dm2 + (w.m2 - dm2) * (omega2 / k2)
        
        d.m2 <- (w.m2 - c.m2 - sm2 - alpha) * m2 +
          a * m1 * (1 - (omega2 / k2))
        
        d.ms2 <- - (c.m2 + alpha) * ms2 + sm2 * m2 +
          a * ms1 * (1 - (omega2 / k2))
        
        d.n1 <- d.f1 + d.m1
        d.ns1 <- d.fs1 + d.ms1
        d.n2 <- d.f2 + d.m2
        d.ns2 <- d.fs2 + d.ms2
        d.N1 <- d.n1 + d.ns1
        d.N2 <- d.n2 + d.ns2
        d.N <- d.N1 + d.N2
        
        list(c(d.f1, d.fs1, d.m1, d.ms1, d.f2, d.fs2, 
               d.m2, d.ms2, d.n1, d.ns1, d.n2, d.ns2,
               d.N1, d.N2, d.N))
      })
    }
    
    init <- c(init['f1'], init['fs1'], 
              init['m1'], init['ms1'],
              init['f2'], init['fs2'], 
              init['m2'], init['ms2'],
              init['n1'], init['ns1'],
              init['n2'], init['ns2'],
              init['N1'], init['N2'], init['N'])
    
    SolveIASA.out <- ode(times = time, 
                         func = SolveIASA.fu, 
                         y = init, 
                         parms = pars,
                         ...)
    
    return(as.data.frame(SolveIASA.out))
  } 
  if (is.null(s.range) & is.null(a.range) & 
        is.null(alpha.range)) {
    output <- SolveIASAfu(pars = pars, init = init, time = time)
    SolveIASA <- list(
      name = 'SolveIASA',
      model = SolveIASAfu,
      pars = pars,
      init = init,
      time = time,
      results = output)
    class(SolveIASA) <- 'capmModels'
    return(SolveIASA)
  } else {
    if(length(a.range) != 2) {
      stop('The length of a.range must be equal to 2.')
    }
    if(length(alpha.range) != 2) {
      stop('The length of alpha.range must be equal to 2.')
    }
    if(length(v.range) != 2) {
      stop('The length of v.range must be equal to 2.')
    }
    if(any(a.range > 1 | alpha.range > 1 | v.range > 1)) {
      stop('Values in a.range, alpha.range,\nv.range and s.range must be lesser or equal to 1.')
    }
    if(any(s.range > 1)) {
      stop('Values in s.range must be lesser or equal to 1.')
    }
    output <- NULL
    paras <- pars
    a.range <- c(a.range[1], pars['a'], a.range[2])
    alpha.range <- c(alpha.range[1], pars['alpha'], alpha.range[2])
    for (i in 1:length(v.range)) {
      for (i1 in 1:length(a.range)) {
        for (i2 in 1:length(alpha.range)) {
          for (i3 in 1:length(s.range)) {
            if (s.fm) {
              paras[c('sf1', 'sm1', 'sf2', 'sm2',
                      'alpha', 'a', 'v')] <- 
                c(s.range[i3], s.range[i3],
                  s.range[i3], s.range[i3],
                  alpha.range[i2], a.range[i1],
                  v.range[i])
            } else {
              paras[c('sf1', 'sf2', 'alpha', 'a', 'v')] <- 
                c(s.range[i3], s.range[i3],
                  alpha.range[i2], a.range[i1], 
                  v.range[i])
            }
            output <- rbind(
              output,
              SolveIASAfu(pars = paras, 
                          init = init, 
                          time = time))
          }
        }
      }
    }
    names(output) <- c(1,2:5, 2:5)
    output <- data.frame(
      rbind(output[, 1:5], output[, c(1, 6:9)]),
      n = c(rowSums(output[, c(2, 4)]), rowSums(output[, c(6, 8)])),
      ns = c(rowSums(output[, c(3, 5)]), rowSums(output[, c(7, 9)])),
      N = c(rowSums(output[, 2:5]), rowSums(output[, 6:9])),
      group = rep(1:2, each = nrow(output)),
      s = rep(s.range, each = length(time)),
      alpha = rep(alpha.range, 
                 each = length(time) * length(s.range)),
      a = rep(a.range, 
                 each = length(time) * length(s.range) * 
                   length(alpha.range)),
      v = rep(v.range, 
               each = length(time) * length(s.range) * 
                 length(alpha.range) * length(a.range)))
    names(output)[1:5] <- c('t', 'f', 'fs', 'm', 'ms')
    SolveIASA <- list(
      name = 'SolveIASA',
      model = SolveIASAfu,
      pars = pars,
      init = init,
      time = time,
      results = output)
    class(SolveIASA) <- 'capmModels'
    return(SolveIASA)
  }
}