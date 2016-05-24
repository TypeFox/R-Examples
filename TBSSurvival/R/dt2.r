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

# t-student density function (an alternative to dt)
.dt2 <- function(x,df) {
  return(exp(lgamma((df+1)/2)-(1/2)*log(df*pi)-lgamma(df/2)-((df+1)/2)*log(1 + (x^2)/df)))
}

# t-student distribution function (an alternative to pt)
.pt2 <- function(x,df) {
  out <- rep(NA,length(x))
  if (!is.null(x[x >= 0]))
    out[x >= 0] <- 1-(1/2)*pbeta(df/(x[x >= 0]^2+df),df/2,1/2)
  if (!is.null(x[x < 0]))
    out[x < 0] <- (1/2)*pbeta(df/(abs(x[x < 0])^2+df),df/2,1/2)
  return(out)
}

# t-student quantile function (an alternative to qt)
.qt2 <- function(q,df) {
  out <- rep(NA,length(q))
  if (!is.null(q[(q >= 0.5) & (q < 1)]))
    out[(q >= 0.5) & (q < 1)] <- sqrt((df/qbeta((-2)*(q[(q >= 0.5) & (q < 1)]-1),df/2,1/2))-df)
  if (!is.null(q[(q < 0.5) & (q > 0)]))
    out[(q < 0.5) & (q > 0)] <- -sqrt((df/qbeta(2*q[(q < 0.5) & (q > 0)],df/2,1/2))-df)
  return(out)
}

# random generate from t-student distribution (an alternative to rt)
.rt2 <- function(n,df) {
  return(.qt2(runif(n),df))
}



