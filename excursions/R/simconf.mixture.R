## simconf.mixture.R
##
##   Copyright (C) 2014 David Bolin, Finn Lindgren
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

simconf.mixture <- function(alpha,
                            mu,
                            Q,
                            w,
                            ind,
                            n.iter=10000,
                            vars,
                            verbose=0,
                            max.threads=0,
                            seed=NULL,
                            mix.samp = TRUE)
{

  if(missing(mu) || !is.list(mu)) {
    stop('Must provide list with mean values')
  } else {
    K <- length(mu)
    n <- length(mu[[1]])
    for(i in seq_len(K)){
      mu[[k]] <- private.as.vector(mu[[k]])
    }
  }
  if(missing(Q) || !is.list(Q)) {
    stop('Must provide list with precision matrices')
  } else {
    if(length(Q) != K){
      stop('Input lists are of different length')
    }
    for(i in seq_len(K)){
      Q[[k]] <- private.as.Matrix(Q[[k]])
    }
  }

  if(missing(w)){
    stop('Must provide list with mixture weights')
  } else {
    w <- private.as.vector(w)
  }
  if(!missing(vars)){
    compute.vars <- FALSE
    if(length(w) != K){
      stop('Input lists are of different length')
    }
    for(i in seq_len(K)){
      vars[[k]] <- private.as.vector(vars[[k]])
    }
  } else {
    compute.vars <- TRUE
    vars <- list()
  }

  if(!missing(ind))
    ind <- private.as.vector(ind)


  if(missing(alpha))
    stop('Must provide significance level alpha')


  if(mix.samp){
    if(missing(ind)){
      ind = seq_len(n)
    }
    Q.chol <- list()
    mu.m <- matrix(0,K,n)
    sd.m <- matrix(0,K,n)
    for(k in seq_len(K))
    {
      Q.chol[[k]] <- chol(Q[[k]])
      mu.m[k,] = mu[[k]]
      if(compute.vars){
        vars[[k]] <- excursions.variances(L = Q.chol[[k]])
      }
      sd.m[k,] = sqrt(vars[[k]])
    }

    limits = c(-1000,1000)
    a.marg = sapply(seq_len(n), function(i) Fmix_inv(p = alpha/2,
                    mu = mu.m[,i], sd = sd.m[,i], w = w, br = limits))

    b.marg = sapply(seq_len(n), function(i) Fmix_inv(p = 1-alpha/2,
                    mu = mu.m[,i], sd = sd.m[,i], w = w, br = limits))

    while(min(a.marg) == limits[1] || max(b.marg) == limits[2])
    {
      limits = 2*limits
      a.marg = sapply(seq_len(n), function(i) Fmix_inv(p = alpha/2,
                      mu = mu.m[,i], sd = sd.m[,i], w = w, br = limits))

      b.marg = sapply(seq_len(n), function(i) Fmix_inv(p = 1-alpha/2,
                      mu = mu.m[,i], sd = sd.m[,i], w = w, br = limits))
    }
    samp <- mix.sample(n.iter,mu,Q.chol,w)
    r.o = optimize(fmix.samp.opt,interval = c(0,alpha),
                   mu=mu.m[,ind], alpha=alpha,
                   sd=sd.m[,ind], w=w, limits = limits,samples=samp[,ind])

    a = sapply(seq_len(n), function(i) Fmix_inv(p = r.o$minimum/2,
               mu = mu.m[,i], sd = sd.m[,i], w = w, br = limits))

    b = sapply(seq_len(n), function(i) Fmix_inv(p = 1-r.o$minimum/2,
               mu = mu.m[,i], sd = sd.m[,i], w = w, br = limits))
  } else {

    if(!missing(ind)){
      if(!is.logical(ind)){
        lind <- rep(FALSE,n)
        lind[ind] = TRUE
        ind <- lind
      }
      cind <- reo <- rep(1,n)
      cind[!ind] = 0
      out <- .C("reordering", nin = as.integer(n), Mp = as.integer(Q[[1]]@p),
                Mi = as.integer(Q[[1]]@i), reo = as.integer(reo),
                cind = as.integer(cind))
      reo = out$reo+1
    } else {
      reo = seq_len(n)
    }
    ireo = NULL
    ireo[reo] = 1:n

    Q.chol <- list()
    mu.m <- matrix(0,K,n)
    sd.m <- matrix(0,K,n)
    for(k in seq_len(K))
    {
      Q.chol[[k]] <- t(as(Cholesky(Q[[k]][reo,reo],perm=FALSE),"Matrix"))
      mu.m[k,] = mu[[k]][reo]
      if(compute.vars){
        vars[[k]] <- excursions.variances(L = Q.chol[[k]])
      }
      sd.m[k,] = sqrt(vars[[k]][reo])
    }

    limits = c(-1000,1000)
    a.marg = sapply(seq_len(n), function(i) Fmix_inv(p = alpha/2,
                    mu = mu.m[,i],sd = sd.m[,i],w = w,br = limits))[ireo]

    b.marg = sapply(seq_len(n), function(i) Fmix_inv(p = 1-alpha/2,
                    mu = mu.m[,i],sd = sd.m[,i],w = w,br = limits))[ireo]

    while(min(a.marg) == limits[1] || max(b.marg) == limits[2])
    {
      limits = 2*limits
      a.marg = sapply(seq_len(n), function(i) Fmix_inv(p = alpha/2,
                      mu = mu.m[,i],sd = sd.m[,i],w = w,br = limits))[ireo]

      b.marg = sapply(seq_len(n), function(i) Fmix_inv(p = 1-alpha/2,
                      mu = mu.m[,i],sd = sd.m[,i],w = w,br = limits))[ireo]
    }

    r.o = optimize(fmix.opt,interval = c(alpha/n,alpha),
                   alpha=alpha,
                   mu=mu.m,
                   sd=sd.m,
                   Q.chol=Q.chol,
                   w=w,
                   ind = ind[reo],
                   limits = limits,
                   max.threads = max.threads,
                   verbose=verbose)

      a = sapply(seq_len(n), function(i) Fmix_inv(p = r.o$minimum/2,
                                                mu = mu.m[,i],
                                                sd = sd.m[,i],
                                                w = w,
                                                br = limits))[ireo]

      b = sapply(seq_len(n), function(i) Fmix_inv(p = 1-r.o$minimum/2,
                                                mu = mu.m[,i],
                                                sd = sd.m[,i],
                                                w = w,
                                                br = limits))[ireo]
  }
  output <- list(a = a[ind],
                 b = b[ind],
                 a.marginal = a.marg[ind],
                 b.marginal = b.marg[ind],
                 mean = mu[ind],
                 vars = vars[ind])
  output$meta = list(calculation="simconf",
                     alpha=alpha,
                     n.iter=n.iter,
                     ind=ind,
                     call = match.call())
  class(output) <- "excurobj"
  return(output)
}

