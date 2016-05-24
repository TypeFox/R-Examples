## excursions.inla.R
##
##   Copyright (C) 2012, 2013, 2014 David Bolin, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

contourmap <- function(mu,
                       Q,
                       vars,
                       n.levels,
                       ind,
                       levels,
                       type = c("standard",
                                "pretty",
                                "equalarea",
                                "P0-optimal",
                                "P1-optimal",
                                "P2-optimal"),
                       compute = list(F=TRUE, measures = NULL),
                       use.marginals=TRUE,
                       alpha,
                       F.limit,
                       n.iter=10000,
                       verbose=FALSE,
                       max.threads=0,
                       seed=NULL)
{
  type <- match.arg(type)

  if(missing(alpha) || is.null(alpha)){
    alpha = 0.1
  }
  if(missing(F.limit)) {
    F.limit = 0.99
  } else {
    F.limit = max(alpha,F.limit)
  }
  if(missing(n.levels) || is.null(n.levels)){
    if(missing(levels) || is.null(levels)){
      stop("Must supply levels or n.levels")
    } else {
      n.levels = length(levels)
    }
  }

  if(!missing(mu))
    mu <- private.as.vector(mu)

  if(!missing(vars))
    vars <- private.as.vector(vars)

  if(!missing(ind))
    ind <- private.as.vector(ind)

  if(!missing(Q))
    Q <- private.as.Matrix(Q)


  measure = NULL
  if(!is.null(compute$measures))
    measure <- match.arg(compute$measures,
                         c("P0", "P1", "P2","P0-bound","P1-bound","P2-bound"),
                         several.ok=TRUE)

  if(type == 'standard')
  {
    if(verbose) cat('Creating contour map\n')
    lp <- excursions.levelplot(mu=mu,n.levels = n.levels,ind = ind,
                               levels = levels,equal.area=FALSE)
  }
  else if(type == 'pretty')
  {
    if(verbose) cat('Creating pretty contour map\n')
    lp <- excursions.levelplot(mu=mu,n.levels = n.levels,ind = ind,
                               levels = levels,equal.area=FALSE,pretty.cm=TRUE)
    n.levels = lp$n.levels
  }
  else if(type == 'equalarea')
  {
    if(verbose) cat('Creating equal area contour map\n')
    lp <- excursions.levelplot(mu = mu,n.levels = n.levels,ind = ind,
                               levels = levels,equal.area=TRUE)
  }
  else if(type == 'P0-optimal' || type == 'P1-optimal' || type == 'P2-optimal')
  {
    if(!missing(levels)){
      warning('Not using supplied levels for optimal contour map\n')
      if(!missing(n.levels)){
        if(n.levels != length(levels)){
          warning('n.levels != length(levels), using n.levels\n')
        }
      } else {
        n.levels = length(levels)
      }
    }
    if(missing(vars) && missing(Q)){
      stop('Variances must be supplied when creating optimal contour map')
    } else if(missing(vars)) {
      vars = excursions.variances(Q=Q)
    }
    if(use.marginals == TRUE){
      if(missing(Q))
        stop('The precision matrix must be supplied unless marginals are used')
    }

    if(type == 'P0-optimal'){
      if(verbose) cat('Creating P0-optimal contour map\n')
      opt.measure = 0
    }else if(type == 'P1-optimal'){
      if(verbose) cat('Creating P1-optimal contour map\n')
      opt.measure = 1
    } else if (type == 'P2-optimal'){
      if(verbose) cat('Creating P2-optimal contour map\n')
      opt.measure = 2
    }

    lp <- excursions.opt.levelplot(mu = mu,vars = vars,Q = Q,
                                   n.levels = n.levels, measure = opt.measure,
                                   use.marginals = use.marginals,ind = ind)
  }

  F.calculated = FALSE
  if(!is.null(measure)){
    if(missing(Q))
      stop('precision matrix must be supplied if measure should be calculated')

    for( i in 1:length(measure)){
      if(measure[i]=="P1") {
        if(n.levels>1){
          if(verbose) cat('Calculating P1-measure\n')
          lp$P1 <- Pmeasure(lp=lp,mu=mu,Q=Q,ind=ind,type=1)
        } else {
          lp$P1 = 1
        }
      } else if(measure[i] == "P2") {
        if(verbose) cat('Calculating P2-measure\n')
        lp$P2 <- Pmeasure(lp=lp,mu=mu,Q=Q,ind=ind,type=2)
      } else if (measure[i] == "P0") {
        if(verbose) cat('Calculating P0-measure and contour map function\n')

        p <- contourfunction(lp=lp, mu=mu,Q=Q ,vars=vars, ind = ind,
                             alpha=alpha, F.limit = F.limit,
                             n.iter=n.iter,max.threads=max.threads,
                             seed=seed,verbose=verbose)
        F.calculated = TRUE
      } else if(measure[i] == "P0-bound"){
        if(missing(vars)){
          vars = excursions.variances(Q=Q)
        }
        lp$P0.bound <- Pmeasure.bound(lp=lp, mu=mu, vars, type=0, ind=ind)
      } else if(measure[i] == "P1-bound"){
        if(missing(vars)){
          vars = excursions.variances(Q=Q)
        }
        lp$P1.bound <- Pmeasure.bound(lp=lp, mu=mu, vars, type=1, ind=ind)
      } else if(measure[i] == "P2-bound"){
        if(missing(vars)){
          vars = excursions.variances(Q=Q)
        }
        lp$P2.bound <- Pmeasure.bound(lp=lp, mu=mu, vars, type=2, ind=ind)
      }
    }
  }
  if(!F.calculated){
    if(is.null(compute$F) || compute$F){
      if(verbose) cat('Calculating contour map function\n')
      p <- contourfunction(lp=lp, mu=mu,Q=Q ,vars=vars, ind = ind,
                           alpha=alpha, F.limit = F.limit,
                           n.iter=n.iter, max.threads=max.threads,
                           seed=seed,verbose=verbose)
      F.calculated = TRUE
    }
  }

  if (missing(ind) || is.null(ind)) {
    ind <- seq_len(length(mu))
  } else if(is.logical(ind)){
    ind <- which(ind)
  }

  if(F.calculated){
    lp$P0 = mean(p$F[ind])
    lp$F = p$F
    lp$E = p$E
    lp$M = p$M
    lp$rho = p$rho
  } else {
    lp$E <- NULL
  }
  lp$meta <- list(calculation="contourmap",
                  F.limit=F.limit,
                  F.computed = compute$F,
                  alpha=alpha,
                  levels=lp$u,
                  type="!=",
                  contourmap.type = type,
                  n.iter=n.iter,
                  mu.range = range(mu[ind]),
                  ind = ind,
                  call = match.call())
  class(lp) <- "excurobj"
  return(lp)
}
