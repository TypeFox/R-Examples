## SA functions
## require(MASS)
SAstep <- function(curr.set, maxvar, 
                   fraction = .3, size.dev = 1) 
{
  if (is.null(curr.set)) { # create a new set
    sample(1:maxvar, round(maxvar * fraction))
  } else { # modify the existing one
    new.size <- length(curr.set)
    if (size.dev > 0) # perhaps change the size
      new.size <- new.size + sample(-size.dev:size.dev, 1)
    if (new.size < 2) new.size <- 2 
    if (new.size > maxvar - 2) new.size <- maxvar - 2

    not.in <- which(!((1:maxvar) %in% curr.set))
    superset <- c(curr.set,
                  sample(not.in, max(size.dev, 1)))
    newset <- sample(superset, new.size)
    ## looks familiar? If yes, then try again
    while (length(newset) == length(curr.set) &&
           !any(is.na(match(newset, curr.set))))
      newset <- sample(superset, new.size)

    newset
  }
}

SAfun <- function(x, response, eval.fun, Tinit, niter = 100,
                  cooling = .05, fraction = .3, ...)
{ # preparations...
  nvar <- ncol(x)
  best <- curr <- SAstep(NULL, nvar, fraction = fraction)
  best.q <- curr.q <- eval.fun(x, response, curr, ...)

  Temp <- Tinit
  for (i in 1:niter) { # Go!
    new <- SAstep(curr, nvar)
    new.q <- eval.fun(x, response, new, ...)
    accept <- TRUE
    if (new.q < curr.q) { # Metropolis criterion
      p.accept <- exp((new.q - curr.q) / Temp)
      if (runif(1) > p.accept) accept <- FALSE
    }
    if (accept) {
      curr <- new
      curr.q <- new.q
      if (curr.q > best.q) { # store best until now
        best <- curr
        best.q <- curr.q
      }
    }
    Temp <- Temp * (1 - cooling)
  }
  list(best = best, best.q = best.q)
}

## the same as SAfun, but now keeping all quality values
SAfun2 <- function(x, response, eval.fun, Tinit, niter = 100,
                   cooling = .05, fraction = .3, ...)
{
  nvar <- ncol(x)
  best <- curr <- SAstep(NULL, nvar, fraction = fraction)
  accepts <- qualities <- rep(NA, niter + 1)
  best.q <- curr.q <- qualities[1] <- eval.fun(x, response, curr, ...)
  accepts[1] <- TRUE
  
  Temp <- Tinit
  for (i in 1:niter) {
    new <- SAstep(curr, nvar)
    new.q <- qualities[i + 1] <- eval.fun(x, response, new, ...)

    accept <- TRUE
    if (new.q < curr.q) {
      p.accept <- exp((new.q - curr.q) / Temp)
      if (runif(1) > p.accept) accept <- FALSE
    }
    accepts[i + 1] <- accept

    if (accept) {
      curr <- new
      curr.q <- new.q
      
      if (curr.q > best.q) {
        best <- curr
        best.q <- curr.q
      }
    }

    Temp <- Temp * (1 - cooling)
  }

  list(best = best, best.q = best.q, qualities = qualities, accepts = accepts)
}

## GA functions

GA.init.pop <- function(popsize, nvar, kmin, kmax)
{
  lapply(1:popsize,
         function(ii, x, min, max) {
           if (min == max) {
             sample(x, min)
           } else {
             sample(x, sample(min:max, 1))
           }},
         nvar, kmin, kmax)
}

## Update July 23: GA.select now also selects variables if none of the
## population members satisfy the min.qlt criterion. Avoids error
## message (reported by Jinsong Zhao)
GA.select <- function(pop, number, qlts, 
                      min.qlt = .4, qlt.exp = 1)
{
  n <- length(pop)
  qlts <- qlts - min(qlts)
  threshold <- quantile(qlts, min.qlt)

  weights <- rep(0.00001, n)
  if(any(OK <- qlts > threshold))
    weights[qlts > threshold] <- 
      qlts[qlts > threshold] ^ qlt.exp

  sample(n, number, replace = TRUE, prob = weights)
}

GA.XO <- function(subset1, subset2)
{
  n1 <- length(subset1)
  n2 <- length(subset2)
  if (n1 == n2) {
    length.out <- n1
  } else {
    length.out <- sample(n1:n2, 1)
  }
  sample(unique(c(subset1, subset2)), length.out)
}

## Update July 23: GA.mut now avoids including more variables than
## maxvar (reported by Jinsong Zhao)
GA.mut <- function(subset, maxvar, mut.prob = .01)
{
  if (runif(1) < mut.prob) { # swap variable in or out
    new <- sample((1:maxvar)[-subset], 1)
    length.out <- min(sample(-1:1 + length(subset), 1),
                      maxvar)# accept at most maxvar variables

    sample(c(new, subset), length.out)
  } else {                   # do nothing
    subset
  }
}

GAfun <- function(X, C, eval.fun, kmin, kmax,
                  popsize = 20, niter = 50,
                  mut.prob = .05, ...)
{
  nvar <- ncol(X) # preparations: the first generation
  pop <- GA.init.pop(popsize, nvar, kmin, kmax)
  pop.q <- sapply(pop,
                  function(subset) eval.fun(X, C, subset, ...))
  best.q <- max(pop.q)
  best <- pop[[which.max(pop.q)]]

  for (i in 1:niter) { # Go!
    new.pop <- 
      lapply(1:popsize, function(j) {
        GA.mut(GA.XO(pop[[GA.select(pop, 1, pop.q)]],
                     pop[[GA.select(pop, 1, pop.q)]]),
               maxvar = nvar, mut.prob = mut.prob)}
             ) # do crossover, and later perhaps mutation
    pop <- new.pop
    pop.q <- sapply(pop,
                    function(subset) eval.fun(X, C, subset, ...))
    if (max(pop.q) > best.q) { # bookkeeping of best solutions
      best.q <- max(pop.q)
      best <- pop[[which.max(pop.q)]]
    }
    if (length(unique(pop.q)) == 1) break # total convergence
  }
  list(best = best, best.q = best.q, n.iter = i)
}

### The same as GAfun, but now keeping best, worst and median quality values

GAfun2 <- function(X, C, eval.fun, kmin, kmax,
                   popsize = 20, niter = 50, mut.prob = .05, ...)
{
  nvar <- ncol(X)
  pop <- GA.init.pop(popsize, nvar, kmin, kmax)
  pop.q <- sapply(pop, function(subset) eval.fun(X, C, subset, ...))
  qualities <- matrix(0, niter + 1, 3)
  qualities[1,] <- quantile(pop.q, c(0, .5, 1))

  best.q <- max(pop.q)
  best <- pop[[which.max(pop.q)]]
    
  for (i in 1:niter) {
    new.pop <-
      lapply(1:popsize, function(j) {
        GA.mut(GA.XO(pop[[GA.select(pop, 1, pop.q)]],
                     pop[[GA.select(pop, 1, pop.q)]]),
               maxvar = nvar, mut.prob = mut.prob)}
             )
    pop <- new.pop
    pop.q <- sapply(pop, function(subset) eval.fun(X, C, subset, ...))
    qualities[i + 1,] <- quantile(pop.q, c(0, .5, 1))
   
    if (max(pop.q) > best.q) {
      best.q <- max(pop.q)
      best <- pop[[which.max(pop.q)]]
    }

    if (length(unique(pop.q)) == 1) break 
  }

  list(best = best, best.q = best.q, n.iter = i,
       qualities = qualities)
}

## evaluation functions
## The dots in lda.loofun are to conform to the format of a generic
## function (which lda.loofun is not...) 
lda.loofun <- function(x, grouping, subset, ...)
{
  lda.obj <- lda(x[, subset, drop = FALSE], 
                 grouping, CV = TRUE)

  sum(lda.obj$class == grouping) - length(subset)
}

pls.cvfun <- function(x, response, subset, ...)
{
  pls.obj <- plsr(response ~ x[,subset],
                  validation = "CV", ...)
  -MSEP(pls.obj, estimate = "CV")$val[pls.obj$ncomp + 1]
}
