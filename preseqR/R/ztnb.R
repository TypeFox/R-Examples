### initial settings of two parameters size and mu in a negative binomial
### distribution for a numeric optimal searching function optim in R
SIZE.INIT = 1
MU.INIT = 0.5

### termination conditions for EM algorithm
TOLERANCE = 1e-10
ITER.TOLERANCE = 1e5


### density function of a truncated zero negative binomial distribution
### size and mu are two parameters for the negative binomial
zerotruncated.dnbinom <- function(x, size, mu, log = FALSE)
{
  ## the density of x in negative binomial
  p <- dnbinom(x, size = size, mu = mu, log = log)

  ## set zeros in x with zero probability
  if (log == FALSE) {
    p[ which(x == 0) ] <- 0
  } else {
    p[ which(x == 0) ] <- -Inf
  }

  ## the density of non-zero in negative binomial
  q <- 1 - dnbinom(0, size = size, mu = mu)

  ## normalize all non-zero values in negrative binomial to generate ZTNB
  if (log == FALSE) {
    return( p/q )
  } else {
    return( p - log(q) )
  }
}


### zerotruncated negative loglikelyhood 
zerotruncated.minus.log.likelyhood <- function(hist.table, size, mu)
{
  prob <- zerotruncated.dnbinom(hist.table[, 1], size, mu, log = TRUE)

  ## negative loglikelyhood
  prob <- -prob
  return( prob %*% hist.table[, 2] )
}


### calculate the negative binomial loglikelyhood
### zero.items is number of items unobserved
### size and mu are parameters in a negative binomial distribution
nb.loglikelyhood <- function(hist.table, zero.items, size, mu)
{
  ## likelyhood of nonzero terms
  log.prob <- dnbinom(hist.table[, 1], size = size, mu = mu, log = TRUE)
  loglikelyhood <- log.prob %*% hist.table[, 2]

  ## add items with zero count
  log.zero.prob <- dnbinom(0, size = size, mu = mu, log = TRUE)
  loglikelyhood <- loglikelyhood + zero.items * log.zero.prob

  return(loglikelyhood)
}


### EM algorithm to fit the histogram with a negative binomial distribution
### hist only includes information for observation
### the number of unobserved items is missing data
preseqR.ztnb.em <- function(n, size=SIZE.INIT, mu=MU.INIT)
{
  checking.hist(n)

  ## setting the number of unobserved items as 0
  zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))

  ## estimate the total number of distinct items
  observed.items <- sum(n[, 2])
  L <- observed.items/( 1 - zero.prob )

  ## expected the number of unobservations
  zero.items <- L*zero.prob

  ## convert zero.items into an integer
  zero.items <- floor(zero.items)

  ## estimated mean and variance
  m <- (n[, 1] %*% n[, 2]) / L
  v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.items )/(L - 1)

  ## target function f
  f <- function(x) {
        return( -nb.loglikelyhood(n, zero.items, size = x, mu = m)/L )
  }

  ## derivative of f
  gr <- function(x) 
  {
    first.term <- ( digamma(x) * zero.items + 
                    digamma(n[, 1] + size) %*% n[, 2] )/L
    second.term <- digamma(x)
    third.term <- log(x) - log(x + m)
    result <- first.term - second.term + third.term
    # f is negative loglikelyhood
    return(-result)
  }

  ## estimate size and mu based on first and second moments
  if (v > m) {
    res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
           lower = 0.0001, upper = 10000)
  } else {
    res <- optim(size, f, gr, method = "L-BFGS-B",
           lower = 0.0001, upper = 10000)
  }

  ## count the times of iteration
  iter <- as.double(1)

  ## initialize the negative loglikelyhood
  loglikelyhood.pre <- Inf

  ## zerotruncated loglikelyhood 
  loglikelyhood <- zerotruncated.minus.log.likelyhood(n, res$par, m)

  ## EM algorithm
  while (( loglikelyhood.pre - loglikelyhood )/observed.items > TOLERANCE &&
           iter < ITER.TOLERANCE)
  {
    ## update negative loglikelyhood
    loglikelyhood.pre <- loglikelyhood

    ## update parameters
    size <- res$par
    mu <- m

### E-step: estimate the number of unobserved items

    ## update the probility an item unobserved
    zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))

    ## estimate the total number of distinct items
    L <- observed.items/( 1 - zero.prob )

    ## update expected number of unobserved items
    zero.items <- L*zero.prob

    ## convert zero.items into an integer
    zero.items <- floor(zero.items)

    ## estimated mean and variance
    m <- (n[, 1] %*% n[, 2])/L
    v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.items )/(L - 1)

### M step: estimate the parameters size and mu
    if (v > m) {
      res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
             lower = 0.0001, upper = 10000)
    } else {
      res <- optim(size, f, gr, method = "L-BFGS-B",
             lower = 0.0001, upper = 10000)
    }
    iter <- iter + 1
    ## zerotruncated loglikelyhood
    loglikelyhood <- zerotruncated.minus.log.likelyhood(n, res$par, m)
  }
  return(list(size = size, mu = mu, loglik = -loglikelyhood.pre))
}


### predict the number of distinct items using EM algorithm
### t is the relative size to inital sample
preseqR.ztnb.estimate <- function(n, t, size=SIZE.INIT, mu=MU.INIT)
{
  checking.hist(n)

  distinct <- sum(n[, 2])

  ## estimate the parameters
  opt <- preseqR.ztnb.em(n, size, mu)
  size <- opt$size
  mu <- opt$mu

  ## the probability of being sampled in the initial experiment
  p <- 1 - dnbinom(0, size = size, mu = mu)

  ## L is the estimated total number of distinct items
  L <- distinct/p

  ## update parameters of negative binomial in the experiment with size n
  mu <- mu*t

  ## the probability of being sampled under the new experiment
  P <- 1 - dnbinom(0, size = size, mu = mu)

  ## return the expected number of new distinct species under the new experiment
  return(L * P - distinct)
}

## predict a complexity curve using EM algorithm
## ss is the step.size
## max.extrapoltion is the maximum value for extrapolation
preseqR.ztnb.species.accum.curve <- function(n, ss = NULL, max.extrapolation = NULL,
                                             size=SIZE.INIT, mu=MU.INIT)
{
  checking.hist(n)

  total.sample <- floor(n[, 1] %*% n[, 2])
  distinct <- sum(n[, 2])

  ## set step.size = total.sample if it is undefined
  if (is.null(ss)) {
    ss <- total.sample
  } else if (ss < 1) {
    write("the step size should be at least one", stderr())
    return(NULL)
  }

  ## set max.extrapolation = 100 * ss if it is undefined
  if (is.null(max.extrapolation)) {

    ## T is the number of experiments; 100 is a magic number
    max.extrapolation <- 100*total.sample
    T <- as.integer( max.extrapolation/ss )

  } else if (max.extrapolation < ss) {
    write("max.extrapolation should be no less then ss", stderr())
    return(NULL)
  } else {
    # T is the number of experiments
    T <- as.integer( max.extrapolation/ss )
  }

  sample.size <- as.double(ss) * (1: T)

  ## estimate parameters
  opt <- preseqR.ztnb.em(n, size, mu)
  size <- opt$size
  mu <- opt$mu

  ## the probability of being sampled in the initial experiment
  p <- 1 - dnbinom(0, size = size, mu = mu)

  ## L is the estimated total number of distinct items
  L <- distinct/p

  ## estimate the item being sampled under new experiments with different size
  t = sample.size/total.sample
  dim(t) = length(t)
  P = apply(t, 1, function(x) {1 - dnbinom(0, size, mu = x * mu)})
  yield.estimates = L*P

  ## combine sample.size and yield.estimates into matrix
  yield.estimates = matrix(c(sample.size, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) = c('sample.size', 'yield.estimate')
  return(yield.estimates)
}


## fitting the negative binoimal distribution to the data by EM algorithm
## ss is the step.size
## max.extrapoltion is the maximum value for extrapolation
## r is a vector of frequencies
preseqR.ztnb.mincount <- function(n, ss = NULL, max.extrapolation = NULL, r=1,
                                  size=SIZE.INIT, mu=MU.INIT)
{
  checking.hist(n)

  total.sample <- floor(n[, 1] %*% n[, 2])
  distinct <- sum(n[, 2])

  ## set step.size = total.sample if it is undefined
  if (is.null(ss)) {
    ss <- total.sample
  } else if (ss < 1) {
    write("the step size should be at least one", stderr())
    return(NULL)
  }

  ## set max.extrapolation = 100 * ss if it is undefined
  if (is.null(max.extrapolation)) {

    ## T is the number of experiments; 100 is a magic number
    max.extrapolation <- 100*total.sample
    T <- as.integer( max.extrapolation/ss )

  } else if (max.extrapolation < ss) {
    write("max.extrapolation should be no less then ss", stderr())
    return(NULL)
  } else {
    # T is the number of experiments
    T <- as.integer( max.extrapolation/ss )
  }

  sample.size <- as.double(ss) * (1: T)

  ## estimate parameters
  opt <- preseqR.ztnb.em(n, size, mu)
  size <- opt$size
  mu <- opt$mu

  ## the probability of being sampled in the initial experiment
  p <- 1 - dnbinom(0, size = size, mu = mu)

  ## L is the estimated total number of distinct items
  L <- distinct/p

  ## estimate the item being sampled under new experiments with different size
  t = sample.size/total.sample
  dim(t) = length(t)
  max.r = max(r)
  yield.estimates = list()
  P = rep(1, length(t))
  for (i in 1:max.r) {
    P = P - apply(t, 1, function(x) dnbinom(i-1, size, mu=x*mu))
    if (i %in% r) {
      yield.estimates = c(yield.estimates, list(L * P))
    }
  }

  #P = lapply(1:max.r, function(x) {
  #      apply(t, 1, function(y) dnbinom(x, size, mu = y * mu))
  #      })
  #yield.estimates = lapply(1:length(r), function(x) {L * P[[x]]})
  yield.estimates = lapply(1:length(r), function(x) {
      estimates <- matrix(c(sample.size, yield.estimates[[x]]), ncol = 2, byrow = FALSE)
      colnames(estimates) <- c("sample.size", paste("yield.estimates(r=", r[x], ")", sep=""))
      estimates
      })

  return(yield.estimates)
}
