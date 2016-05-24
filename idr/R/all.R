# Author: Qunhua Li
# Affiliation: UC Berkeley, Peter Bickel's group
# Reference:
#   Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)
#   Measuring reproducibility of high-throughput experiments.
#   Annals of Applied Statistics. In press.
#
# Upon using this code, you automatically agree to cite the reference above
# in all publications that include results generated from this code or this
# R code itself.
#
# License: GPL-2 

est.IDR <-
function(x, mu, sigma, rho, p, eps=0.001, max.ite=30){

# This is the main function estimates parameters from the copula mixture models
# Input:
#   x: a m x n numeric matrix, where m is the number of replicates and n is
#      the number of entries on each replicate
#   mu, sigma, rho, p: starting values for mean, sd, correlation coefficient
#      and mixing proportion of the reproducible component
#
# Output:
#   (p, rho, mu, sigma): estimates
#   loglik: the mixture likelihood at the end of iterations
#   bic: bic
#   e.z: proterior probability for each observation to be reproducible, i.e.
#        1-local idr
#   loglik.trace: the trace of mixture likelihood  
#
# Implementation details: iterate between the following two steps:
# 1. raw values are first transformed into pseudovalues
# 2. EM is used to compute the underlining structure, which is a mixture
#    of two normals
  
    conv <- function(old,new) abs(new-old) < eps*(1+abs(new))
    
    x1 <- x[, 1]
    x2 <- x[, 2]
    x1.cdf.func <- ecdf(x1)
    x2.cdf.func <- ecdf(x2)
    afactor <- length(x1)/(length(x1) + 1)
    x1.cdf <- x1.cdf.func(x1) * afactor
    x2.cdf <- x2.cdf.func(x2) * afactor
    para <- list()
    para$mu <- mu
    para$sigma <- sigma
    para$rho <- rho
    para$p <- p
    j <- 1
    to.run <- TRUE
    loglik.trace <- c()
    loglik.inner.trace <- c()
    z.1 <- get.pseudo.mix(x1.cdf, para$mu, para$sigma, para$rho, 
        para$p)
    z.2 <- get.pseudo.mix(x2.cdf, para$mu, para$sigma, para$rho, 
        para$p)
    while (to.run) {
        i <- 1
        while (to.run) {
            e.z <- e.step.2normal(z.1, z.2, para$mu, para$sigma, 
                para$rho, para$p)
            para <- m.step.2normal(z.1, z.2, e.z)
            if (i > 1) 
                l.old <- l.new
            l.new <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, 
                para$rho, para$p)
            loglik.inner.trace[i] <- l.new
            if (i > 1) {
                to.run <- !conv(loglik.inner.trace[i-1],loglik.inner.trace[i])
            }
            i <- i + 1
        }
       
        z.1 <- get.pseudo.mix(x1.cdf, para$mu, para$sigma, para$rho, 
            para$p)
        z.2 <- get.pseudo.mix(x2.cdf, para$mu, para$sigma, para$rho, 
            para$p)
        if (j > 1) 
            l.old.outer <- l.new.outer
        l.new.outer <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, 
            para$rho, para$p)

        loglik.trace[j] <- l.new.outer
        if (j == 1) 
            to.run <- TRUE
        else {
            if (j > max.ite) 
                to.run <- FALSE
            else to.run <- !conv(l.old.outer,l.new.outer)
            
        }
        j <- j + 1
    }
    idr <- 1 - e.z
    o <- order(idr)
    idr.o <- idr[o]
    idr.rank <- rank(idr.o, ties.method = "max")
    top.mean <- function(index, x) {
        mean(x[1:index])
    }
    IDR.o <- sapply(idr.rank, top.mean, idr.o)
    IDR <- idr
    IDR[o] <- IDR.o
    return(list(para = list(p = para$p, rho = para$rho, mu = para$mu, 
        sigma = para$sigma), loglik = l.new, loglik.trace = loglik.trace, 
        idr = 1 - e.z, IDR = IDR))
}

#est.IDR.hist <- function(x, mu, sigma, rho, p, eps, max.ite=30){

  # currently it is a placeholder, to add real function
  
#  est <- est.IDR(x, mu, sigma, rho, p, eps, max.ite=30)
#  return(est)
#}


#est.IDR.discrete <- function(x, mu, sigma, rho, p, eps, max.ite=30){

  # currently it is a placeholder, to add real function
#  est <- est.IDR(x, mu, sigma, rho, p, eps, max.ite=30)
#  return(est)
#}


select.IDR <-
  function(x, IDR.x, IDR.level){

    # Select observations that exceeding a given IDR level
    #
    # Input:
    #   x: a m by n numeric matrix, where m= num of replicates,
    #      n=num of observations. Numerical values representing the
    #      significance of the observations, where signals are expected to have
    #      large values, for example, -log(p-value).   Currently, m=2. 
    #   IDR.x: Irreproducibile discovery rate for each entry of x.
    #      It is computed from est.IDR().
    #   IDR.level: IDR cutoff, a numerical value between [0,1].
    #
    # Output:
    #   x: Observations that are selected.
    #   n: Number of observations that are selected.
    #   IDR.level: IDR cutoff, a numerical value between [0,1].
    #

    is.selected <- IDR.x < IDR.level
    x.selected <- x[is.selected,]
    n.selected <- nrow(x.selected)
    
    return(list(x=x.selected, n=n.selected, IDR.level=IDR.level))
  }


get.correspondence <-
function(x1, x2, t, spline.df=NULL){
  
# compute the correspondence profile
# Input:
#   x1: Data values or ranks of the data values on list 1, a vector of
#    numeric values. Large values need to be significant signals. If small
#    values represent significant signals, rank the signals reversely
#    (e.g. by ranking negative values) and use the rank as x1.
#   x2: Data values or ranks of the data values on list 2, a vector of
#    numeric values. Large values need to be significant signals. If small
#    values represent significant signals, rank the signals reversely
#    (e.g. by ranking negative values) and use the rank as x1.
#   t: A numeric vector between 0 and 1 in ascending order. t is the
#    right-tail percentage.
#   spline.df: Degree of freedom for spline, to control the smoothness of the
#    smoothed curve. 
#
# Output:
#  psi: the correspondence profile in terms of the scale of percentage,
#     i.e. between (0, 1)
#  dpsi: the derivative of the correspondence profile in terms of the scale
#     of percentage, i.e. between (0, 1)}
#  psi.n: the correspondence profile in terms of the scale of the number of
#     observations
#  dpsi.n: the derivative of the correspondence profile in terms of the scale
#     of the number of observations}
#
#  Each object above is a list consisting of the following items:
#  t: upper percentage (for psi and dpsi) or number of top ranked
#  observations (for psi.n and dpsi.n) 
#  value: psi or dpsi
#  smoothed.line: smoothing spline
#  ntotal: the number of observations
#  jump.point: the index of the vector of t such that psi(t[jump.point])
#              jumps up due to ties at the low values. This only
#	      happends when data consists of a large number of discrete
#	      values, e.g. values imputed for observations
#	      appearing on only one replicate. 
  
  psi.dpsi <- get.uri.2d(x1, x2, t, t, spline.df)

  psi <- list(t=psi.dpsi$tv[,1], value=psi.dpsi$uri, 
              smoothed.line=psi.dpsi$uri.spl, ntotal=psi.dpsi$ntotal, 
              jump.point=psi.dpsi$jump.left)

  dpsi <- list(t=psi.dpsi$t.binned, value=psi.dpsi$uri.slope,
               smoothed.line=psi.dpsi$uri.der, ntotal=psi.dpsi$ntotal, 
               jump.point=psi.dpsi$jump.left)

  psi.n <- list(t=psi$t*psi$ntotal, value=psi$value*psi$ntotal, 
                 smoothed.line=list(x=psi$smoothed.line$x*psi$ntotal, 
                                    y=psi$smoothed.line$y*psi$ntotal),
                 ntotal=psi$ntotal, jump.point=psi$jump.point) 

  dpsi.n <- list(t=dpsi$t*dpsi$ntotal, value=dpsi$value, 
                 smoothed.line=list(x=dpsi$smoothed.line$x*dpsi$ntotal, 
                                    y=dpsi$smoothed.line$y),
                 ntotal=dpsi$ntotal, jump.point=dpsi$jump.point) 

  return(list(psi=psi, dpsi=dpsi, psi.n=psi.n, dpsi.n=dpsi.n))
}


get.uri.2d <-
function(x1, x2, tt, vv, spline.df=NULL){

  o <- order(x1, x2, decreasing=TRUE)
  
  # sort x2 by the order of x1
  x2.ordered <- x2[o]
  
  tv <- cbind(tt, vv)
  ntotal <- length(x1) # number of peaks    

  uri <- apply(tv, 1, comp.uri, x=x2.ordered)

  # compute the derivative of URI vs t using small bins
  uri.binned <- uri[seq(1, length(uri), by=4)]
  tt.binned <- tt[seq(1, length(uri), by=4)]
  uri.slope <- (uri.binned[2:(length(uri.binned))] - uri.binned[1:(length(uri.binned)-1)])/(tt.binned[2:(length(uri.binned))] - tt.binned[1:(length(tt.binned)-1)])

  # smooth uri using spline
  # first find where the jump is and don't fit the jump
  # this is the index on the left
  # jump.left.old  <- which.max(uri[-1]-uri[-length(uri)])
  short.list.length <- min(sum(x1>0)/length(x1), sum(x2>0)/length(x2))

  if(short.list.length < max(tt)){
    jump.left <- which(tt>short.list.length)[1]-1
  } else {
    jump.left <- which.max(tt)
  }

#  reversed.index <- seq(length(tt), 1, by=-1)
#  nequal <- sum(uri[reversed.index]== tt[reversed.index])
#  temp  <- which(uri[reversed.index]== tt[reversed.index])[nequal]
#  jump.left <- length(tt)-temp
 
  if(jump.left < 6){
   jump.left <- length(tt)
  }
    
 
  if(is.null(spline.df))
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=6.4)
  else{
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=spline.df)
  }
  # predict the first derivative
  uri.der <- predict(uri.spl, tt[1:jump.left], deriv=1)

  invisible(list(tv=tv, uri=uri, 
                 uri.slope=uri.slope, t.binned=tt.binned[2:length(uri.binned)], 
                 uri.spl=uri.spl, uri.der=uri.der, jump.left=jump.left,
                 ntotal=ntotal))
 }



#plot.correspondence <-
#  function(corr.profile, plot.type, plot.all){

   # Plot the correspondence profile for the output from get.correspondence()
   #  
   # Input:
   #   corr.profile: correspondence profile computed from get.correspondence()
   #   plot.type: "uri", "duri", "uri.n", "duri.n"
   #     "uri": upper rank intersection (y-axis) vs t (x-axis), where t is
   #            the percentage from upper side. This is the psi plot in the
   #            reference paper. 
   #     "duri": derivative of upper rank intersaction (y-axis) vs t (x-axis).
   #            This is the psi' plot in the reference paper
   #     "uri.n": same as "uri", except that x-axis is the number of
   #            observations
   #     "duri.n": same as "duri", except that x-axis is the number of
   #            observations
   #   plot.all: logical value.  
   #     If the rank lists consist of a large number of ties at the bottom
   #    (e.g. the same low value is imputed to the list for the observations
   #    that appear on only one list), it may be desirable to plot only
   #    observations before hitting the ties. True: plot all observations.
   #    False: plot observations before hitting the ties. 
    
#  if(plot.type=="uri"){  

    # plot correspondence curve on the scale of percentage
#    if(plot.all){
      
#      plot(corr.profile$psi$t, corr.profile$psi$value, xlab="t", ylab="psi", xlim=c(0, max(corr.profile$psi$t)), ylim=c(0, max(corr.profile$psi$value)), cex.lab=2)
#      lines(corr.profile$psi$smoothed.line, lwd=4)
#      abline(coef=c(0,1), lty=3)
#    } else {
    
    # If the rank lists consist of a large number of ties at the bottom
    # (e.g. the same low value is imputed to the list for the observations
    # that appear on only one list), it may be desirable to plot only
    # observations before hitting the ties. Then it can be plotted using the
    # following
#      plot(corr.profile$psi$t[1:corr.profile$psi$jump.point], corr.profile$psi$value[1:corr.profile$psi$jump.point], xlab="t", ylab="psi", xlim=c(0, max(corr.profile$psi$t[1:corr.profile$psi$jump.point])), ylim=c(0, max(corr.profile$psi$value[1:corr.profile$psi$jump.point])), cex.lab=2)
#      lines(corr.profile$psi$smoothed.line, lwd=4)
#      abline(coef=c(0,1), lty=3)
#    }
    
#  } else {

#    if(plot.type=="duri"){
      
      # plot the derivative of correspondence curve on the scale of percentage 
#      plot(corr.profile$dpsi$t, corr.profile$dpsi$value, xlab="t", ylab="psi'", xlim=c(0, max(corr.profile$dpsi$t)), ylim=c(0, max(corr.profile$dpsi$value)), cex.lab=2)
#      lines(corr.profile$dpsi$smoothed.line, lwd=4)
#      abline(h=1, lty=3)
#    } else {

#      if(plot.type=="uri.n"){
        
        # plot correspondence curve on the scale of the number of observations
#        plot(corr.profile$psi.n$t, corr.profile$psi.n$value, xlab="t", ylab="psi", xlim=c(0, max(corr.profile$psi.n$t)), ylim=c(0, max(corr.profile$psi.n$value)), cex.lab=2)
#        lines(corr.profile$psi.n$smoothed.line, lwd=4)
#        abline(coef=c(0,1), lty=3)
#      } else {
        
        # plot the derivative of correspondence curve on the scale of the
        # number of observations
#        plot(corr.profile$dpsi.n$t, corr.profile$dpsi.n$value, xlab="t", ylab="psi'", xlim=c(0, max(corr.profile$dpsi.n$t)), ylim=c(0, max(corr.profile$dpsi.n$value)), cex.lab=2)
#        lines(corr.profile$dpsi.n$smoothed.line, lwd=4)
#        abline(h=1, lty=3)
#      }
#    }
#  }
  
#}


##################### start internal functions #############

comp.uri <-
function(tv, x){

  # Internal function
  #
  # Compute the value of Psi and Psi' for a given t
  # An internal function to compute Psi for a given (t, v)
  
  # Input:
  #   tv: A vector of two numeric values, t and v. Both t and v are in [0,1]
  #   x: A numeric vector of x, sorted by the order of y
  #
  # Output:
  #   A numeric value of Psi(t, v)

  n <- length(x)
  qt <- quantile(x, prob=1-tv[1]) # tv[1] is t
  sum(x[1:ceiling(n*tv[2])] >= qt)/n

}


d.binormal <-
  function(z.1, z.2, mu, sigma, rho){

  # Internal function
  #  
  # Compute the log-density for parameterized bivariate Gaussian
  # distribution N(mu, mu, sigma, sigma, rho).
  #
  # Input:
  #   z.1: a numerical data vector on coordinate 1. 
  #   z.2: a numerical data vector on coordinate 1.
  #   mu: mean 
  #   sigma: standard deviation 
  #   rho: correlation coefficient
  #  
  # Output:
  #   Log density of bivariate Gaussian distribution
  #    N(mu, mu, sigma, sigma, rho).
  

  loglik <- (-log(2)-log(pi)-2*log(sigma) - log(1-rho^2)/2 - (0.5/(1-rho^2)/sigma^2)*((z.1-mu)^2 -2*rho*(z.1-mu)*(z.2-mu) + (z.2-mu)^2))

  return(loglik)
}

e.step.2normal <-
function(z.1, z.2, mu, sigma, rho, p){

  # Internal function
  #
  # Expectation step in the EM algorithm for parameterized bivariate
  # 2-component Gaussian mixture models with (1-p)N(0, 0, 1, 1, 0) +
  # pN(mu, mu, sigma, sigma, rho).

  # Input:
  #  z.1: a numerical data vector on coordinate 1.
  #  z.2: a numerical data vector on coordinate 2.
  #  mu: mean for the reproducible component.
  #  sigma: standard deviation of the reproducible component.
  #  rho: correlation coefficient of the reproducible component.
  #  p: mixing proportion of the reproducible component.
  #
  # Output:
  #  e.z: a numeric vector, where each entry represents the estimated expected
  #  conditional probability that an observation is in the reproducible
  #  component.

  e.z <- p/((1-p)*exp(d.binormal(z.1, z.2, 0, 1, 0)-d.binormal(z.1, z.2, mu, sigma, rho))+ p)
  
  invisible(e.z)
}

get.pseudo.mix <-
function(x, mu, sigma, rho, p){

  # Internal function
  #
  # Compute the pseudo values of a mixture model from the empirical CDF
  #
  # Input:
  #   x: A vector of values of empirical CDF
  #   mu: Mean of the reproducible component in the mixture model on the
  #       latent space
  #   sigma: Standard deviation of the reproducible component in the mixture
  #       model on the latent space
  #   rho: Correlation coefficient of the reproducible component in the
  #       mixture model on the latent space
  #   p: Mixing proportion of the reproducible component in the mixture model
  #      on the latent space
  #
  # Output:
  #   The values of a mixture model corresponding to the empirical CDF
  
  # first compute cdf for points on the grid
  # generate 200 points between [-3, mu+3*sigma]
  nw <- 1000
  w <- seq(min(-3, mu-3*sigma), max(mu+3*sigma, 3), length=nw) 
  w.cdf <- p*pnorm(w, mean=mu, sd=sigma) + (1-p)*pnorm(w, mean=0, sd=1)

  i <- 1

  quan.x <- rep(NA, length(x))

  for(i in c(1:nw)){
    index <- which(x >= w.cdf[i] & x < w.cdf[i+1])
    quan.x[index] <- (x[index]-w.cdf[i])*(w[i+1]-w[i])/(w.cdf[i+1]-w.cdf[i]) +w[i]
  }

  index <- which(x < w.cdf[1])
  if(length(index)>0)
    quan.x[index] <- w[1]

  index <- which(x > w.cdf[nw])
  if(length(index)>0)
    quan.x[index] <- w[nw]  
  
  invisible(quan.x)
}

loglik.2binormal <-
function(z.1, z.2, mu, sigma, rho, p){

  # Internal function
  #
  # Compute the log-likelihood for parameterized bivariate 2-component Gaussian
  # mixture models with (1-p)N(0, 0, 1, 1, 0) + pN(mu, mu, sigma,sigma, rho).
  # 
  # Input:
  #  z.1: a numerical data vector on coordinate 1. 
  #  z.2: a numerical data vector on coordinate 1. 
  #  mu: mean for the reproducible component. 
  #  sigma: standard deviation of the reproducible component. 
  #  rho: correlation coefficient of the reproducible component. 
  #  p: mixing proportion of the reproducible component. 

  # Output:
  #  Log-likelihood of the bivariate 2-component Gaussian mixture models
  #  $(1-p)N(0, 0, 1, 1, 0) + N(mu, mu, sigma, sigma, rho)$.
  
  l.m <- sum(d.binormal(z.1, z.2, 0, 1, 0)+log(p*exp(d.binormal(z.1, z.2, mu, sigma, rho)-d.binormal(z.1, z.2, 0, 1, 0))+(1-p)))
  
  return(l.m) 
}


m.step.2normal <-
function(z.1, z.2, e.z){

  # Internal function
  # Maximization step in the EM algorithm for parameterized bivariate
  # 2-component Gaussian mixture models with $(1-p)N(0, 0, 1, 1, 0) +
  # pN(mu, mu, sigma^2, sigma^2, rho)$.
  #
  # Input:
  #  z.1: a numerical data vector on coordinate 1.
  #  z.2: a numerical data vector on coordinate 2.
  #  e.z: a vector of expected conditional probability that the $i$th
  #  observation is reproducible. 
  #
  # Output:
  #  Estimated parameters, basically a list including elements 
  #  p: the estimated mixing proportion of the reproducible component.
  #  mu: the estimated mean for the reproducible component.
  #  sigma: the estimated standard deviation of the reproducible component.
  #  rho: the estimated correlation coefficient of the reproducible component.
  
  
  p <- mean(e.z)
  mu <- sum((z.1+z.2)*e.z)/2/sum(e.z) 
  sigma <- sqrt(sum(e.z*((z.1-mu)^2+(z.2-mu)^2))/2/sum(e.z))
  rho <- 2*sum(e.z*(z.1-mu)*(z.2-mu))/(sum(e.z*((z.1-mu)^2+(z.2-mu)^2)))

  return(list(p=p, mu=mu, sigma=sigma, rho=rho))
}

select.IDR <-
  function(x, IDR.x, IDR.level){

    is.selected <- IDR.x < IDR.level
    x.selected <- x[is.selected,]
    n.selected <- nrow(x.selected)
    
    return(list(x=x.selected, n=n.selected, IDR.level=IDR.level))
  }

