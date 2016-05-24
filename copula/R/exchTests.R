## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' Test of exchangeability for bivariate EV copulas based on the Pickands
##' or CFG estimators -- see SJS paper
##'
##' @title Test of exchangeability based on An
##' @param x the data
##' @param N number of multiplier replications
##' @param estimator Pickands or CFG
##' @param derivatives based on "An" or "Cn"
##' @param m grid size
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
exchEVTest <- function(x, N = 1000, estimator = "CFG", derivatives = "Cn", m = 100)
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- pobs(x)

  ## make grid
  g <- seq(1/m, 0.5, len = m)

  ## compute the test statistic
  s <- .C(evsymtest_stat,
          as.double(-log(u[,1])),
          as.double(-log(u[,2])),
          as.integer(n),
          as.double(g),
          as.integer(m),
          as.integer(estimator == "CFG"),
          stat = double(1))$stat

  if (derivatives == "Cn")
    s0 <- .C(evsymtest,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(estimator == "CFG"),
             as.integer(N),
             s0 = double(N))$s0
  else
    s0 <- .C(evsymtest_derA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(estimator == "CFG"),
             as.integer(N),
             s0 = double(N))$s0

  structure(class = "htest",
	    list(method = paste("Test of exchangeability for bivariate extreme-value copulas with argument 'estimator' set to '",
                 estimator,"', argument 'derivatives' set to '",derivatives,"' and argument 'm' set to ", m, sep=""),
                 statistic = c(statistic = s),
                 p.value = (sum(s0 >= s) + 0.5) / (N + 1),
                 data.name = deparse(substitute(x))))
}

##' Test of exchangeability for bivariate copulas based on the
##' empirical copula -- see SJS paper
##'
##' @title Test of exchangeability based on Cn
##' @param x the data
##' @param N number of multiplier replications
##' @param m grid size; if 0, use pseudo-observations
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
exchTest <- function(x, N = 1000, m = 0)
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- pobs(x)

  ## make grid
  if (m > 0)
    {
      xis <- yis <- seq(1/m, 1 - 1/m, len = m)
      g <- as.matrix(expand.grid(xis, yis))
      ng <- m^2
    }
  else
    {
      g <- u
      ng <- n
    }

  ## compute the test statistic
  s <- .C(exchtestCn_stat,
          as.double(u[,1]),
          as.double(u[,2]),
          as.integer(n),
          as.double(g[,1]),
          as.double(g[,2]),
          as.integer(ng),
          stat = double(1))$stat

  s0 <- .C(exchtestCn,
           as.double(u[,1]),
           as.double(u[,2]),
           as.integer(n),
           as.double(g[,1]),
           as.double(g[,2]),
           as.integer(ng),
           as.integer(N),
           s0 = double(N))$s0

  structure(class = "htest",
	    list(method = paste("Test of exchangeability for bivariate copulas with argument 'm' set to",m),
                 statistic = c(statistic = s),
                 p.value = (sum(s0 >= s) + 0.5) / (N + 1),
                 data.name = deparse(substitute(x))))
}
