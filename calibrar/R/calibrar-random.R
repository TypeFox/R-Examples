# Random numbers ----------------------------------------------------------

# main function for generate random numbers, semi-optimize for speed
rtnorm2 = function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf, lq=0.025) {
  
  if(length(mean)!=length(sd)) stop("mean and sd must have the same length.")
  if(length(lower)!=length(upper)) stop("upper and lower must have the same length.")
  if(length(mean)!=length(upper)) stop("mean/sd and lower/upper must have the same length.")
  
  if(all(lower==-Inf) & all(upper==Inf)) {
    return(.rnorm(n=n, mean=mean, sd=sd))
  }
  
  output = matrix(nrow=length(mean), ncol=n)
  
  ind = (lower == -Inf & upper == +Inf)

  if(any(ind)) {
    output[ind, ] = .rnorm(n=n, mean=mean[ind], sd=sd[ind])
    output[!ind,] = .rtnorm2(n=n, mean=mean[!ind], sd=sd[!ind], 
                             lower=lower[!ind], upper=upper[!ind])
  } else {
    output = .rtnorm2(n=n, mean=mean, sd=sd, lower=lower, upper=upper)
  }
  
  return(output)
}

# auxiliar functions for random numbers generation
.rnorm = function(n, mean=0, sd=1) {
  
  if(length(mean)!=length(sd)) stop("mean and sd must have the same length.")
  output = matrix(rnorm(n*length(mean), mean=mean, sd=sd), nrow=length(mean))
  return(output)
}

.rnorm2 = function(n, mean=0, sd=1, lower=-Inf, upper=Inf) {
  
  if(length(mean)!=length(sd)) stop("mean and sd must have the same length.")
  if(length(lower)!=length(upper)) stop("upper and lower must have the same length.")
  if(length(mean)!=length(upper)) stop("mean/sd and lower/upper must have the same length.")
  
  output = matrix(rnorm(n*length(mean), mean=mean, sd=sd), nrow=length(mean))
  
  ind1 = which(output < lower, arr.ind=TRUE)[, 1]
  if(length(ind1)>0) {
    ind2 = which(output < lower)
    output[ind2] = lower[ind1]
  }
  
  ind1 = which(output > upper, arr.ind=TRUE)[, 1]
  if(length(ind1)>0) {
    ind2 = which(output > upper)
    output[ind2] = upper[ind1]
  }
  
  return(output)
  
}

.rtnorm = function(n, mean=0, sd=1, lower = -Inf, upper = Inf) {
  
  if(length(mean)!=length(sd)) stop("mean and sd must have the same length.")
  if(length(lower)!=length(upper)) stop("upper and lower must have the same length.")
  if(length(mean)!=length(upper)) stop("mean/sd and lower/upper must have the same length.")
  
  output = matrix(.rtnormx(n*length(mean), mean=mean, sd=sd, lower=lower, upper=upper), 
                  nrow=length(mean))
  return(output)
}

.rtnorm2 = function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf, lq=0.025) {
  
  if(length(mean)!=length(sd)) stop("mean and sd must have the same length.")
  if(length(lower)!=length(upper)) stop("upper and lower must have the same length.")
  if(length(mean)!=length(upper)) stop("mean/sd and lower/upper must have the same length.")
  
  ind = pnorm(q=lower, mean=mean, sd=sd) - pnorm(q=upper, mean=mean, sd=sd) + 1 > lq
  
  if(all(!ind)) {
    output = .rnorm2(n=n, mean=mean, sd=sd, lower=lower, upper=upper)
    return(output)
  }
  
  if(all(ind)) {
    output = .rtnorm(n=n, mean=mean, sd=sd, lower=lower, upper=upper)
    return(output)
  }
  
  output = matrix(nrow=length(mean), ncol=n)
  output[!ind, ] = .rnorm2(n=n, mean=mean[!ind], sd=sd[!ind], 
                           lower=lower[!ind], upper=upper[!ind])
  output[ind, ]  = .rtnorm(n=n, mean=mean[ind], sd=sd[ind], 
                           lower=lower[ind], upper=upper[ind])
  
  return(output)
  
}

# borrowed from msm package, not optimized function, to be updated
.rtnormx = function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  
  if (length(n) > 1) n <- length(n)
  mean <- rep(mean, length = n)
  sd <- rep(sd, length = n)
  lower <- rep(lower, length = n)
  upper <- rep(upper, length = n)
  lower <- (lower - mean)/sd
  upper <- (upper - mean)/sd
  ind <- seq(length = n)
  ret <- numeric(n)
  alg <- ifelse(lower > upper, -1, 
                ifelse(((lower < 0 & upper == Inf) | 
                          (lower == -Inf & upper > 0) | (is.finite(lower) & 
                           is.finite(upper) & (lower < 0) & (upper > 0) & (upper - lower > sqrt(2 * pi)))), 0, 
                       ifelse((lower >= 0 & 
                                 (upper > lower + 2 * sqrt(exp(1))/(lower + sqrt(lower^2 + 4))*
                                    exp((lower * 2 - lower * sqrt(lower^2 + 4))/4))),1, 
                              ifelse(upper <= 0 & (-lower > -upper + 2 * sqrt(exp(1))/(-upper + sqrt(upper^2 + 4)) * 
                                                     exp((upper * 2 - -upper * sqrt(upper^2 +4))/4)), 2, 3))))
  
  ind.nan <- ind[alg == -1]
  ind.no <- ind[alg == 0]
  ind.expl <- ind[alg == 1]
  ind.expu <- ind[alg == 2]
  ind.u <- ind[alg == 3]
  
  ret[ind.nan] <- NaN
  
  while (length(ind.no) > 0) {
    y <- rnorm(length(ind.no))
    done <- which(y >= lower[ind.no] & y <= upper[ind.no])
    ret[ind.no[done]] <- y[done]
    ind.no <- setdiff(ind.no, ind.no[done])
  }
  stopifnot(length(ind.no) == 0)
  
  while (length(ind.expl) > 0) {
    a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4))/2
    z <- rexp(length(ind.expl), a) + lower[ind.expl]
    u <- runif(length(ind.expl))
    done <- which((u <= exp(-(z - a)^2/2)) & (z <= upper[ind.expl]))
    ret[ind.expl[done]] <- z[done]
    ind.expl <- setdiff(ind.expl, ind.expl[done])
  }
  stopifnot(length(ind.expl) == 0)
  
  while (length(ind.expu) > 0) {
    a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4))/2
    z <- rexp(length(ind.expu), a) - upper[ind.expu]
    u <- runif(length(ind.expu))
    done <- which((u <= exp(-(z - a)^2/2)) & (z <= -lower[ind.expu]))
    ret[ind.expu[done]] <- -z[done]
    ind.expu <- setdiff(ind.expu, ind.expu[done])
  }
  stopifnot(length(ind.expu) == 0)
  
  while (length(ind.u) > 0) {
    z <- runif(length(ind.u), lower[ind.u], upper[ind.u])
    rho <- ifelse(lower[ind.u] > 0, exp((lower[ind.u]^2 - 
                                           z^2)/2), ifelse(upper[ind.u] < 0, exp((upper[ind.u]^2 - 
                                                                                    z^2)/2), exp(-z^2/2)))
    u <- runif(length(ind.u))
    done <- which(u <= rho)
    ret[ind.u[done]] <- z[done]
    ind.u <- setdiff(ind.u, ind.u[done])
  }
  stopifnot(length(ind.u) == 0)
  
  return(ret * sd + mean)
  
}
