# Source Functions

# ## plotting utilities
# # inputs not necessary to be sorted; they are sorted internally
# my.radial <- function(r, theta, ...){
#   radial.plot(c(r[order(theta)]), c(theta[order(theta)]),
#               rp.type = "p", show.grid.label = TRUE, radial.lim = c(0, 0.5),
#               ...)
# }
# # rotate a matrix
# rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees
# # # reference: http://stackoverflow.com/questions/16496210/rotate-a-matrix-in-r

# general ellipse function to generate the true boundary
# (a, b) - semidiameter paramters (a > b)
# (r0, theta0) are the relative polar coodinates for the center to the origin (center of the image)
# phi - the rotation angle of the a-axis relative to the polar axis
ellipse <- function(a, b, r0 = 0, theta0 = 0, phi = 0){
  function(theta){
    P = r0 * ((b^2 - a^2) * cos(theta + theta0 - 2 * phi) + (a^2 + b^2) * cos(theta - theta0))
    R = (b^2 - a^2) * cos(2 * theta - 2 * phi) + a^2 + b^2
    Q = sqrt(2) * a * b * sqrt(R - 2 * r0^2 * sin(theta - theta0)^2)
    r = (P + Q)/R
    return(r)
  }
}

# regular triangle
# arg - h is the height (altitude) of the regular triangle

triangle <- function(h){

  # theta is from 0 to 2*pi
  # for scalar input
  triangle.scalar <- function(theta, h){

    if (theta >= pi/3 & theta < 2*pi/3){
      theta <- -theta + 2 * pi/3
    }

    if (theta >= 2 * pi/3 & theta < 4 * pi/3){
      theta <- theta - 2 * pi/3
    }

    if (theta >= 4 * pi/3 &  theta <= 2 * pi){
      theta <- 2 * pi - theta
    }

    ret = (2 * h / 3) /(cos(theta) + sqrt(3) * sin(theta))
    return(ret)

  }

  ret = function(theta) c(sapply(theta, function(theta) triangle.scalar(theta, h)))
  return(ret)
}



# from parameters to observation; set seed before usage
# m: number of pixels each direction: total number n = m^2
# pi.in: In probability to be 1;
# pi.out: Out probability to be 1;
# design = 'D' (deterministic) or 'J'(Jitteredly Random Design), or 'C' (Completely Random Design)
# gamma.fun - equation of the function 'gamma' in a polar system, i.e. a function of theta
# polar coordinate: (theta, radius)
par2obs <- function(m, pi.in, pi.out, design, gamma.fun){
  # center of obs is (0.5, 0.5)
  center = c(0.5, 0.5)

  obs <- matrix(NA, m, m) # pre-locate observation
  if (design == 'D'){
    x.axis = (col(obs) - 1)/m + 1/(2*m)
    y.axis = (m - row(obs))/m + 1/(2*m)
  }

  if (design == "J"){
    x.axis = (col(obs) - 1)/m + 1/(2*m) + runif(m^2, min = -1/(2 * m), max = 1/(2 * m))
    y.axis = (m - row(obs))/m + 1/(2*m) + runif(m^2, min = -1/(2 * m), max = 1/(2 * m))
  }

  # if (design == "C"){
  #   x.axis = matrix(runif(m^2), m, m)
  #   y.axis = matrix(runif(m^2), m, m)
  # }
  # need to be sorted if "C" is used; so far, we don't use it

  r.obs = sqrt((x.axis - center[1])^2 + (y.axis - center[2])^2)
  theta.obs <- atan2(y.axis - 1/2, x.axis - 1/2) #[-pi, pi]
  # atan2(y, x) returns the angle between (x, y) but has the range [-pi, pi]
  # we want to have the range [0, 2*pi) - note: exclude 2*pi since it's the same as 0
  theta.obs[theta.obs < 0 ] = theta.obs[theta.obs < 0 ] + 2*pi  # [0, 2*pi)
  # we won't see c(0, 2 * pi) in 'theta.obs' since the center c(1/2, 1/2) doesn't share the x-axis or y-axis with any observation locations

  obsLabel = (r.obs < gamma.fun(theta.obs)) # label - inside or outside
  n.In = sum(obsLabel)
  n.Out = sum(!obsLabel)
  obs[obsLabel] = rbinom(n.In, size = 1, prob = pi.in)
  obs[!obsLabel] = rbinom(n.Out, size = 1, prob = pi.out)

  return(list(intensity = obs, theta.obs = theta.obs, r.obs = r.obs,
              center = center, x = x.axis, y = y.axis))
  # not sure it's necessary to save all those; how about also save (x.axis, y.axis)?
  # come back here to find out what to save later on
  # so far, we can forget about the "pixel" world; now we have (location, observation)
  # the pixel is used when plot an image - but even it is not necessary to use "pixel" even for this purpose
  # we can add point iteratively by adding (location, observation)
}


# eigenfunctions: sqrt(2*pi), sqrt(pi)*cosx, sqrt(pi)*sinx, ...
# the 1st, 2nd, 3rd, ... eigenfunctions
eigen.fun = function(n){

  k1 = n %% 2;
  k2 = (n - k1)/2;

  if (n == 1) {
    ret = function(x) 1/sqrt(2*pi) + 0 * x
  }

  if (n > 1){
    if (k1 == 0){
      ret = function(x) 1/sqrt(pi) * cos(k2 * x)
    }

    if (k1 == 1){
      ret = function(x) 1/sqrt(pi) * sin(k2 * x)
    }
  }

  ret
}



## Univariate slice sampling by Neal (2003)
# UNIVARIATE SLICE SAMPLING WITH STEPPING OUT AND SHRINKAGE.
#
# Performs a slice sampling update from an initial point to a new point that
# leaves invariant the distribution with the specified log density function.
#
# Arguments:
#
#   x0    Initial point
#   g     Function returning the log of the probability density (plus constant)
#   w     Size of the steps for creating interval (default 1)
#   m     Limit on steps (default infinite)
#   lower Lower bound on support of the distribution (default -Inf)
#   upper Upper bound on support of the distribution (default +Inf)
#   gx0   Value of g(x0), if known (default is not known)
#


uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{

#   uni.slice.calls = 0
#   uni.slice.evals = 0
#   # Keep track of the number of calls made to this function.
#
#   uni.slice.calls <- uni.slice.calls + 1

  # Find the log density at the initial point, if not already known.

  if (is.null(gx0))
  { # uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }

  # Determine the slice level, in log terms.

  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      # uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      # uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      # uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      # uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }

  # Shrink interval to lower and upper bounds.

  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  {
    x1 <- runif(1,L,R)

    # uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)

    if (gx1>=logy) break

    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.

  attr(x1,"log.density") <- gx1
  return (x1)

}

BayesBD.binary <- function(obs, ini.mean = 0.4, n.run = 10000, n.burn = 1000,
                           J = 10, output.all = FALSE){

  ### MCMC setup ----
  # use smoothed means as the prior mean & initial estimate: the constant 'ini.mean'
  mu.fun <- function(theta)  theta * 0 + ini.mean
  mu = c(mu.fun(obs$theta.obs))

  # number of eigenvalues and eigenfunctions
  L = 2 * J + 1 # include the intercept

  # hyper parameters: (alpha.tau, beta.tau) and (alpha.lambda, beta.lambda) are for (shape, rate)
  alpha.tau = 500; beta.tau = 1; alpha.lambda = 2; beta.lambda = 1


  # global calculation
  obs$demean = obs$r.obs - mu
  obs$idxOne = which(obs$intensity == 1) # only applicable for binary images
  tmp.mat = sapply(1:L, function(k) eigen.fun(k)(c(obs$theta.obs)))

  ## evaluate the matrix at theta.plot
  #tmp.mat.theta.plot = sapply(1:L, function(k) eigen.fun(k)(c(theta.plot)))
  ## evaluate initial mu function at theta.plot
  #mu.theta.plot = mu.fun(theta.plot)


  # record samples
  an.smp = matrix(NA, L, n.run)
  pi.smp = matrix(NA, 2, n.run) # 1st row - pi.in
  tau.smp = rep(NA, n.run)
  lambda.smp = rep(NA, n.run)

  # initial values
  an.ini = rep(0, L)
  pi.in.ini = sum((obs$r.obs < mu) * obs$intensity)/sum(obs$r.obs < mu)
  pi.out.ini = sum((obs$r.obs > mu) * obs$intensity)/sum(obs$r.obs > mu)
  tau.ini = 500
  lambda.ini = 1

  # related computation of the initial
  # all the eigenvalues at lambda.ini
  cn.ini = sapply(0:J, function(ith) 2*pi*exp(-2 * lambda.ini) * besselI(2*lambda.ini, ith))
  eigen.cn.ini = rep(cn.ini, each = 2)[-1]

  # used in the likelihood
  alpha.ini = log(pi.in.ini * (1 - pi.out.ini) / (pi.out.ini * ( 1 - pi.in.ini)))
  beta.ini = log((1 - pi.in.ini)/(1 - pi.out.ini))


  # r - r.hat
  diff.ini = obs$demean - c(tmp.mat %*% an.ini) # initialize
  tmp = (diff.ini < 0)
  n.in <- sum(tmp)
  n.in.one <- sum(tmp[obs$idxOne]) # more efficient

  ## MCMC update

  # return the conditional posterior log density of an[k] given an[-k] and the rest of pars
  # diff = observed radius - estimated radius
  # important to have 'diff' here to save computation!
  # use 'attr()' to save the estimation of the log density
  logL.an.k <- function(k, an, alpha, beta, eigen.cn, tau, diff, tmp.mat){
    eigen.cn.k = eigen.cn[k]
    an.k = an[k]
    tmp.mat.k = c(tmp.mat[,k])

    ret = function(an.k.candidate){
      diff.new = diff - (an.k.candidate - an.k) * tmp.mat.k
      tmp = (diff.new <= 0)
      n.in <- sum(tmp)
      # n.in.one <- sum(tmp * obs$intensity) #require obs
      n.in.one <- sum(tmp[obs$idxOne])
      ret2 = n.in.one * alpha + n.in * beta - an.k.candidate^2 / eigen.cn.k * tau/2
      attr(ret2, 'diff') <- diff.new # not sure whether this is a burden or not -  to check
      attr(ret2, 'n.in') <- n.in
      attr(ret2, 'n.in.one') <- n.in.one
      return(ret2)
    }

    return(ret)
  }

  logL.lambda = function(tau, an){
    ret = function(lambda){
      cn = sapply(0:J, function(ith) 2*pi*exp(-2 * lambda) * besselI(2*lambda, ith))
      eigen.cn = rep(cn, each = 2)[-1]
      ret2 = -1/2 * (sum(log(eigen.cn)) + tau * sum(an^2/eigen.cn)) + (alpha.lambda - 1) * log(lambda) - beta.lambda * lambda
      return(ret2)
    }
    return(ret)
  }

  start <- proc.time()

  # inside the loop:
  # update alpha.ini, beta.ini, diff.ini, n.in, n.in.one for each ith to save computation
  for (ith in 1:n.run){
    # update an one by one
    for (k in 1:L){
      f.tmp = logL.an.k(k, an.ini, alpha.ini, beta.ini,
                        eigen.cn.ini, tau.ini, diff.ini, tmp.mat)

      tmp2delete <- uni.slice(an.ini[k], g = f.tmp,
                              gx0 = n.in.one * alpha.ini + n.in * beta.ini - an.ini[k]^2 / eigen.cn.ini[k] * tau.ini/2)

      # diff.ini.t1 = diff.ini - tmp.mat[,k] * (tmp2delete - an.ini[k])
      tmp2delete.2 = attr(tmp2delete, 'log.density')
      diff.ini = attr(tmp2delete.2, 'diff')
      n.in = attr(tmp2delete.2, 'n.in')
      n.in.one = attr(tmp2delete.2, 'n.in.one')
      an.ini[k] <- tmp2delete
    }


    # update tau
    a.star = alpha.tau + L/2
    b.star = beta.tau + sum(an.ini^2 /eigen.cn.ini)/2
    tau.ini = rgamma(1, shape = a.star, rate = b.star)

    # update pi.in, pi.out
    r.mat = c(tmp.mat %*% an.ini + mu)
    tmp <- (obs$r.obs <= r.mat) # require obs.info
    n.in <- sum(tmp)
    n.in.one <- sum(tmp * obs$intensity) #require obs
    n.out = length(r.mat) - n.in
    n.out.one = sum(obs$intensity) - n.in.one
    x1 = rbeta(1, n.in.one, n.in - n.in.one)
    x2 = rbeta(1, n.out.one, n.out - n.out.one)
    pi.in.ini = max(x1, x2)
    pi.out.ini = min(x1, x2)

    # update alpha.ini and beta.ini
    # used in the likelihood
    alpha.ini = log(pi.in.ini * (1 - pi.out.ini) / (pi.out.ini * ( 1 - pi.in.ini)))
    beta.ini = log((1 - pi.in.ini)/(1 - pi.out.ini))


    # update lambda
    lambda.ini = uni.slice(lambda.ini, g = logL.lambda(tau = tau.ini, an = an.ini), lower = 0, upper = 500)

    # save samples
    an.smp[,ith] = an.ini
    pi.smp[,ith] = c(pi.in.ini, pi.out.ini)
    lambda.smp[ith] = lambda.ini
    tau.smp[ith] = tau.ini

    # report
    if (ith %% floor(n.run/100) == 1){
      duration = proc.time() - start;
      cat(ith, "th iteration:", "ETA:", duration[3]/ith * (n.run - ith), "\n")
    }
  }

  Total.Time = proc.time() - start

  # Estimates
  an.MC = apply(an.smp[, n.burn:n.run], 1, mean)

  gamma.MC <- function(theta.arg){
    c(sapply(1:L, function(k) eigen.fun(k)(theta.arg)) %*% an.MC + mu.fun(theta.arg))
  }

  # return
  retsult = list(gamma.hat = gamma.MC)

  if (output.all){
    result = list(an.smp = an.smp, pi.smp = pi.smp,
                  tau.smp = tau.smp, lambda.smp = lambda.smp,
                  gamma.hat = gamma.MC)
  }

  return(result)
}


