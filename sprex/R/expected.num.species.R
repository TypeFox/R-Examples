#' @name expected.num.species
#' @title Expected Number of Species
#' @description Calculate the expected number of species for a given 
#'   sample size.
#' 
#' @param m number of samples.
#' @param f a vector of species frequencies where \code{f[i]} is the number 
#'   of species represented by only \code{i} samples.
#' @param f0.func a function that computes the number of unobserved 
#'   species (f0).
#' @param ... other arguments to \code{f0.func}.
#' 
#' @return a vector or matrix (depending on whether \code{m} is a scalar or vector, 
#'   respectively) of the estimated number of species (\code{s.ind}) seen 
#'   in \code{m} samples, and the standard deviation (\code{sd.s.ind}).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references Eqns 4, 5, 9, and 10 in Colwell, R.K., A. Chao, N.J. Gotelli, 
#'   S.-Y. Lin, C.X. Mao, R.L. Chazdon, and J.T. Longino. 2012. Models and 
#'   estimators linking individual-based and sample-based rarefaction, 
#'   extrapolation and comparison of assemblages. 
#'   Journal of Plant Ecology 5(1):3-21.
#'  
#' @examples
#' data(osa.old.growth)
#' f <- expand.freqs(osa.old.growth)
#' expected.num.species(60, f = f, f0.func = Chao1)
#' 
#' expected.num.species(c(60, 70, 75), f = f, f0.func = Chao1)
NULL

#' @export
#' 
.s.ind.n.m <- function(f0, f1, n, m.star, s.obs) {
  if(f0 == 0) return(s.obs)
  # calculate Sind(n + m*) Eqn 9
  # term.1 <- f1 / (n * f0)
  # term.2 <- (1 - term.1) ^ m.star
  term.2 <- exp(-(m.star / n) * (f1 / f0))
  s.obs + (f0 * (1 - term.2))
}

#' @rdname expected.num.species
#' @export
#' 
expected.num.species <- function(m, f, f0.func, ...) {
  x <- f0.func(f, ...)
  s.est <- unname(x["s.est"])
  s.obs <- unname(x["s.obs"])
  f0 <- unname(x["f0"])
  n <- unname(x["n"])

  result <- sapply(m, function(m.i) {
    s.ind <- if(m.i <= n) {
      # calculate Sind(m) Eqn 4 & 5
      a.km <- sapply(1:length(f), function(k) {
        # have to simplify binomial coefficient and use log-transforms
        # to avoid NaNs for large values of n
        if(k <= (n - m.i)) {
          num <- (n - k - m.i + 1):(n - m.i)
          denom <- (n - k + 1):n
          exp(sum(log(num)) - sum(log(denom)))
        } else 0
      })
      s.ind <- s.obs - sum(a.km * f)
      term.1 <- ((1 - a.km) ^ 2) * f
      term.2 <- s.ind ^ 2 / s.est
      var.s.ind <- sum(term.1) - term.2
      if(var.s.ind < 0) var.s.ind <- 0
      c(s.ind = s.ind, sd.s.ind = sqrt(var.s.ind))
    } else if(f0 == 0) {
      # if f0 = 0, Eqn 9 produces NaNs
      c(s.ind = s.obs, sd.s.ind = 0) 
    } else {
      m.star <- m.i - n
      s.ind <- .s.ind.n.m(f0, f[1], n, m.star, s.obs)
      
      # variance estimation using partial derivatives from 
      #   Colwell et al 2012 Eqn 10 (modified from code by Alex Curtis)
      calc.dS.df <- function(i, f) {
        f.m1 <- f.p1 <- f
        f.m1[i] <- max(c(0, f.m1[i] - 1))
        f.p1[i] <- f.p1[i] + 1
        f0.m1 <- f0.func(f.m1, ...)
        f0.p1 <- f0.func(f.p1, ...)
        dS.m1 = .s.ind.n.m(f0.m1["f0"], f.m1[1], n, m.star, f0.m1["s.obs"])
        dS.p1 = .s.ind.n.m(f0.p1["f0"], f.p1[1], n, m.star, f0.p1["s.obs"])
        (dS.p1 - dS.m1) / 2
      }
      
      dS.df <- rep(NA, length(f))
      dS.df[1] <- calc.dS.df(1, f)
      if(length(dS.df) > 1) {
        dS.df[2] <- calc.dS.df(2, f)
        if(length(dS.df) > 2) {
          dS.df.gt.2 <- calc.dS.df(3, f)
          for(i in 3:length(dS.df)) dS.df[i] <- dS.df.gt.2
        }
      }
      
      var.s.ind <- if(length(f) > 1) {
        var.off.diag <- sum(sapply(1:(length(f) - 1), function(i) {
          sum(sapply((i + 1):length(f), function(j) {
            dS.df[i] * dS.df[j] * (-f[i] * f[j] / s.est)
          }))
        }))
        var.diag <- sum(sapply(1:length(f), function(i) {
          dS.df[i] ^ 2 * (f[i] * (1 - f[i] / s.est))
        }))
        var.diag + (2 * var.off.diag)
      } else {
        dS.df[1] ^ 2 * (f[1] * (1 - f[1] / s.est))
      }
      
      c(s.ind = s.ind, sd.s.ind = sqrt(var.s.ind))
    }
    
    c(s.ind, m = m.i, x)
  })
  
  drop(t(result))
}
