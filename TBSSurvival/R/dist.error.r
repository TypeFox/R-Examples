# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2012-2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
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

#######################################################################
## builder of the four functions for well-known distributions
## d(x,xi), p(x,xi), q(x,xi), r(x,xi) for the chosen distribution and
## the last element in the list is the distribution name.
dist.error <- function(dist="norm") {
  if ((dist != "norm") && (dist != "doubexp") && (dist != "cauchy") &&
      (dist != "t")    && (dist != "logistic") && (dist != "logis2") && (dist != "all"))
    stop("TBS: Distribution not available at dist.error")

  switch(dist,
         ## Normal distribution
         norm = list(
           d = function(x,xi) dnorm(x,mean=0,sd=sqrt(xi)), # density
           p = function(x,xi) pnorm(x,mean=0,sd=sqrt(xi)), # distr
           q = function(x,xi) qnorm(x,mean=0,sd=sqrt(xi)), # quantile
           r = function(x,xi) rnorm(x,mean=0,sd=sqrt(xi)), # generation
           name = "norm"
           ),
         ## t-student distribution
         t = list(
           d = function(x,xi) .dt2(x,df=xi), # density
           p = function(x,xi) .pt2(x,df=xi), # distr
           q = function(x,xi) .qt2(x,df=xi), # quantile
           r = function(x,xi) .rt2(x,df=xi), # generation
           name = "t"
           ),
         ## Cauchy distribution
         cauchy = list(
           d = function(x,xi) dcauchy(x,location=0,scale=xi), # density
           p = function(x,xi) pcauchy(x,location=0,scale=xi), # distr
           q = function(x,xi) qcauchy(x,location=0,scale=xi), # quantile
           r = function(x,xi) rcauchy(x,location=0,scale=xi), # generation
           name = "cauchy"
           ),
         ## Laplace/Double exponential distribution
         doubexp  = list(
           d = function(x,xi) dnormp(x,sigmap=xi,mu=0,p=1), # density
           p = function(x,xi) pnormp(x,sigmap=xi,mu=0,p=1), # distr
           q = function(x,xi) qnormp(x,sigmap=xi,mu=0,p=1), # quantile
           r = function(x,xi) rnormp(x,sigmap=xi,mu=0,p=1), # generation
           name = "doubexp"
           ),
         ## Logistic distribution
         logis2 = list(
           d = function(x,xi) .dlogis2(x,s=xi), # density
           p = function(x,xi) .plogis2(x,s=xi), # distr
           q = function(x,xi) .qlogis2(x,s=xi), # quantile
           r = function(x,xi) .rlogis2(x,s=xi), # generation
           name = "logis2"
           ),
         logistic = list(
           d = function(x,xi) dlogis(x,location=0,scale=xi), # density
           p = function(x,xi) plogis(x,location=0,scale=xi), # distr
           q = function(x,xi) qlogis(x,location=0,scale=xi), # quantile
           r = function(x,xi) rlogis(x,location=0,scale=xi), # generation
           name = "logistic"
           ),
         all = list(dist.error("norm"),dist.error("t"),dist.error("cauchy"),dist.error("doubexp"),dist.error("logistic"))
         )
}
