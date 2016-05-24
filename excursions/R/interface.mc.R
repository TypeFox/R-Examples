## interface.mc.R
##
##   Copyright (C) 2015 David Bolin
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

contourmap.mc <- function(samples,
                       n.levels,
                       ind,
                       levels,
                       type = c("standard",
                                "equalarea",
                                "P0-optimal",
                                "P1-optimal",
                                "P2-optimal"),
                       compute = list(F=TRUE, measures = NULL),
                       alpha,
                       verbose=FALSE)
{

  if(missing(samples)){
    stop("Must supply samples.")
  } else {
    samples <- as(samples,"matrix")
  }
  if(!missing(ind))
    ind <- private.as.vector(ind)

  mu <- rowMeans(samples)
  type <- match.arg(type)

  if(missing(alpha) || is.null(alpha)){
    alpha = 0.1
  }

  measure = NULL
  if(!is.null(compute$measures))
    measure <- match.arg(compute$measures,
                         c("P0", "P1", "P2"),
                         several.ok=TRUE)

  if(type == 'standard')
  {
    if(verbose) cat('Creating contour map\n')
    lp <- excursions.levelplot(mu=mu,n.levels = n.levels,ind = ind,
                               levels = levels,equal.area=FALSE)
  }
  else if(type == 'equalarea')
  {
    if(verbose) cat('Creating equal area contour map\n')
    lp <- excursions.levelplot(mu = mu,n.levels = n.levels,ind = ind,
                               levels = levels,equal.area=TRUE)
  }
  else if(type == 'P0-optimal' || type == 'P1-optimal' || type == 'P2-optimal')
  {
    warning('Pk-optimal contour maps not implemented, using standard.\n')
    lp <- excursions.levelplot(mu=mu,n.levels = n.levels,ind = ind,
                               levels = levels,equal.area=FALSE)

  }

  F.calculated = FALSE
  if(!is.null(measure)){

    for( i in 1:length(measure)){
      if(measure[i]=="P1") {
          if(verbose) cat('Calculating P1-measure\n')
          lp$P1 <- Pmeasure.mc(lp=lp,mu=mu,X=samples,ind=ind,type=1)
      } else if(measure[i] == "P2") {
        if(verbose) cat('Calculating P2-measure\n')
        lp$P2 <- Pmeasure.mc(lp=lp,mu=mu,X=samples,ind=ind,type=2)
      } else if (measure[i] == "P0") {
        if(verbose) cat('Calculating P0-measure and contour map function\n')

        p <- contourfunction.mc(lp=lp, mu=mu,X=samples, ind = ind,
                             alpha=alpha, verbose=verbose)
        F.calculated = TRUE
      }
    }
  }
  if(!F.calculated){
    if(is.null(compute$F) || compute$F){
      if(verbose) cat('Calculating contour map function\n')
      p <- contourfunction.mc(lp=lp, mu=mu,X=samples, ind = ind,
                           alpha=alpha, verbose=verbose)
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
                  F.limit=0,
                  alpha=alpha,
                  levels=lp$u,
                  type="!=",
                  n.iter=dim(samples)[2],
                  mu.range = range(mu[ind]),
                  ind = ind)
  class(lp) <- "excurobj"
  return(lp)
}

simconf.mc <- function(samples,
                       alpha,
                       ind,
                       verbose=FALSE)
{

  if(missing(samples))
    stop('Must provide matrix with samples')

  if(missing(alpha))
    stop('Must provide significance level alpha')

  if(missing(ind)){
    ind = seq_len(dim(samples)[1])
  }

  a.marg = apply(samples,1,quantile,1,probs=c(alpha/2))
  b.marg = apply(samples,1,quantile,1,probs=c(1-alpha/2))

  #Simple golden section search
  lb = 0
  ub = alpha
  gr = 2/(sqrt(5) + 1)
  x1 = ub - gr*(ub - lb)
  x2 = lb + gr*(ub - lb)


  f1 = fsamp.opt(x1,samples=samples[ind,], verbose=verbose)
  f2 = fsamp.opt(x2,samples=samples[ind,], verbose=verbose)

  while (abs(ub - lb) > 1e-4) {
    if (f2 < 1-alpha) {
      # optimum is to the left of x2
      ub = x2
      x2 = x1
      f2 = f1
      x1 = ub- gr*(ub- lb)
      f1 = fsamp.opt(x1,samples=samples[ind,], verbose=verbose)
    } else {
      lb = x1
      x1 = x2
      f1 = f2
      x2 = lb + gr*(ub - lb)
      f2 = fsamp.opt(x2,samples=samples[ind,], verbose=verbose)
    }
  }

  rho  = (lb + ub)/2
  cat(rho)
  a = apply(samples,1,quantile,1,probs=c(rho/2))
  b = apply(samples,1,quantile,1,probs=c(1-rho/2))

  return(list(a = a[ind],
              b = b[ind],
              a.marginal = a.marg[ind],
              b.marginal = b.marg[ind]))

}



excursions.mc <- function(samples,
                          alpha,
                          u,
                          type,
                          rho,
                          reo,
                          ind,
                          max.size,
                          verbose=FALSE)
{

  if(missing(alpha))
    stop('Must specify error probability')

  if(missing(u))
    stop('Must specify level')

  mu = rowMeans(samples)

  if(missing(type))
    stop('Must specify type of excursion set')

  if(!missing(ind) && !missing(reo))
    stop('Either provide a reordering using the reo argument or provied a set of nodes using the ind argument, both cannot be provided')


  F.limit = 1

  if(verbose)
    cat("Calculate marginals\n")
  marg <- excursions.marginals.mc(X=samples, type = type, rho = rho,
                                  mu = mu, u = u)

  if (missing(max.size)){
    m.size = length(mu)
  } else {
    m.size = max.size
  }
  if (!missing(ind)) {
    if(is.logical(ind)){
      indices = ind
      if(missing(max.size)){
        m.size = sum(ind)
      } else {
        m.size = min(sum(ind),m.size)
      }
    } else {
      indices = rep(FALSE,length(mu))
      indices[ind] = TRUE
      if(missing(max.size)){
        m.size = length(ind)
      } else {
        m.size = min(length(ind),m.size)
      }
    }
  } else {
    indices = rep(TRUE,length(mu))
  }

  if(verbose)
    cat("Calculate permutation\n")
  if(missing(reo)){
    reo <- excursions.permutation(marg$rho, indices, use.camd = FALSE)
  }

  if(verbose)
    cat("Calculate limits\n")
  limits <- excursions.setlimits(marg=marg, type=type,u=u,
                                 mu=rep(0,length(mu)),QC=FALSE)


  res = mcint(X=samples[reo,],a=limits$a[reo],b=limits$b[reo])

  n = length(mu)
  ii = which(res$Pv[1:n] > 0)
  if (length(ii) == 0) i=n+1 else i=min(ii)

  F = Fe  = E = G = rep(0,n)
  F[reo] = res$Pv
  Fe[reo] = res$Ev
  ireo = NULL
  ireo[reo] = 1:n

  ind.lowF = F < 1-F.limit
  E[F>1-alpha] = 1

  if(type == '=') {
    F=1-F
  }

  if(type == "<") {
    G[mu>u] = 1
  } else {
    G[mu>=u] = 1
  }

  F[ind.lowF] = Fe[ind.lowF] = NA

  M = rep(-1,n)
  if (type=="<") {
    M[E==1] = 0
  } else if (type == ">") {
    M[E==1] = 1
  } else if (type == "!=" || type == "=") {
    M[E==1 & mu>u] = 1
    M[E==1 & mu<u] = 0
  }

  if (missing(ind) || is.null(ind)) {
    ind <- seq_len(n)
  } else if(is.logical(ind)){
    ind <- which(ind)
  }
  vars <- rowSums((samples-rowMeans(samples))^2)/(dim(samples)[2]-1)
  output <- list(F = F,
                 G = G,
                 M = M,
                 E = E,
                 mean = mu,
                 vars=vars,
                 rho=marg$rho,
                 meta=(list(calculation="excursions",
                            type=type,
                            level=u,
                            F.limit=F.limit,
                            alpha=alpha,
                            n.iter=dim(samples)[2],
                            method='MC',
                            ind=ind,
                            reo=reo,
                            ireo=ireo,
                            Fe=Fe,
                            LDL=FALSE)))
  class(output) <- "excurobj"
  output
}

