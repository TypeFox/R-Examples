# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# logistic density function (an alternative to dlogis)
.dlogis2 <- function(x,s) {
  return(exp(-(x/s))/(s*(1+exp(-x/s))^2))
}

# logistic distribution function (an alternative to plogis)
.plogis2 <- function(x,s) {
  return(1/(1+exp(-x/s)))
}

# logistic quantile function (an alternative to qlogis)
.qlogis2 <- function(q,s) {
  return(-s*log(1/q-1))
}

# random generate from logistic distribution (an alternative to rlogis)
.rlogis2 <- function(n,s) {
  return(.qlogis2(runif(n),s))
}



