rt.star.k <- function(psi.unique.values, m.star, m.theta, dprime)
 {
  ## First compute some quantities like the location and scale
  ## parameters of the t-distribution.
  aa.prime <- dprime[1]; bb.prime <- dprime[2]
  cc.prime <- dprime[3]; dd.prime <- dprime[4]
  half.m.theta <- 0.5 * m.theta
  t.location <- cc.prime
  t.scale <- sqrt(bb.prime * dd.prime / aa.prime)
  t.df <- 2.0 * aa.prime
  ## print(paste("t.location", t.location, "t.scale", t.scale, "t.df", t.df))
  ##
  ## Now compute the t cdf at psi.unique.
  ## Assure psi.unique is sorted!
  psi.unique <- sort(psi.unique.values)
  ## print(paste("psi.unique", psi.unique))
  t.cdf <- pt((psi.unique - t.location) / t.scale, t.df)
  ## psi.unique will partition the real line into m.star + 1 disjoint
  ## intervals. Compute the probabilities assigned by the t distribution
  ## to each of the intervals.
  log.t.probs <- log(c(t.cdf, 1.0) - c(0.0, t.cdf))
  ## print(paste("log.t.probs", log.t.probs))
  ##
  ## Next compute the value of the K function in the m.star + 1 intervals
  ## A bit of mathematics shows that C * K where C is a constant is more
  ## tractable than K itself for numerical accuracy. C can be chosen to be
  ## 1/Gamma(half.m.theta)^2 in which case things simplify considerably.
  ## Then we populate the K values by exploiting its symmetry.
  k.values <- 0:m.star
  offset <- 2 * lgamma(half.m.theta)
  for (i in 0:m.star) {
    k.values[i+1] <- lgamma(half.m.theta + i) + lgamma(half.m.theta + m.star - i) - offset
  }
  ## print(paste("k.values", k.values))
  ##
  ## Now again for numerical accuracy, we center the log of the
  ## probabilities
  log.product <- log.t.probs - k.values
  ## print(paste("log.product", log.product))
  centered.log.product <- log.product - mean(log.product)
  ## print(paste("centered.log.product", centered.log.product))
  probs <- exp(centered.log.product)
  ## print(paste("rt probs", probs))
  index <- sample(0:m.star, 1, prob = probs)
  if (index == 0) {
    b <- t.cdf[1]
    a <- 0.0
  } else if (index == m.star) {
    b <- 1.0
    a <- t.cdf[m.star]
  } else {
    b <- t.cdf[index + 1]
    a <- t.cdf[index]
  }
  t.scale  * qt( a + (b - a) * runif(n=1), t.df ) + t.location
}

### The Gibbs sampler using dirichlet mixtures prior
dirichlet.c <- function(data, ncycles=10, M=1,d=c(.1,.1,0,1000),
                        start=NULL)
  {
  cl <- match.call()
  inv.gam.par <- d[1:2]; norm.par <- d[3:4]
  if (!is.numeric(data) || !is.matrix(data))
    stop("data must be numeric matrix")
  psi.hat <- data[,1]; se.psi.hat <- data[,2]
  if (!all(se.psi.hat >= 0)) {
     stop("negative standard error in column 2 of data.\n")
  }
  if (M <= 0)  {
   stop("Dirichlet precision parameter M must be greater than zero.\n")
  }
  if (!is.numeric(d) || length(d) != 4){
    stop("hyperparameter must be numeric of length 4")
  }
  else if (d[1] <= 0 || d[2] <= 0 || d[4] <= 0)
  {
     stop("Gamma shape, Gamma scale, and normal variance hyperparameters
          must all be greater than zero\n")
  }
  names(d) <- c("gam.shape","gam.scale","mu","var.factor")
  aa <- d[1]; bb <- d[2]; cc <- d[3]; dd <- d[4]  
  cycle.number <- 0
  nstudies <- length(psi.hat)    
  if (is.null(start)) {
    start.user <- FALSE
    psi.vec.current <- psi.hat
    psi.current <- mean(psi.vec.current)
    tau.current <- sqrt(var(psi.vec.current))
    start <- c(psi.vec.current,psi.current,tau.current)
  }
  else if (!is.null(start)){
    if(!is.numeric(start) || !(length(start)==nstudies+2)) {
     stop("initial param must be num. vec. length nstudies+2, or omit")
    }
    else {
     start.user <- TRUE
     psi.vec.current <- start[1:nstudies]
     psi.current <- start[nstudies+1]
     tau.current <- start[nstudies+2]
    }
  }  
  vec.to.save  <- c(psi.vec.current, psi.current, tau.current)
  output.matrix <- matrix(0, ncycles+1, nstudies+2)
  output.matrix[1,] <- vec.to.save
  study.names <- dimnames(data)[[1]]
  if (!is.null(study.names))
    {
    dimnames(output.matrix) <- list(NULL,param=
                                  c(study.names,"mu","tau"))
    }
  else
    {dimnames(output.matrix) <- list(cycle.id=0:ncycles,param=
                               c(as.character(1:nstudies),"mu","tau"))
    }
  for (ic in (1:ncycles)) {
                                        #    print(paste("IC = ", ic))
    for (i in (1:nstudies)) {
                                        #      print(paste("I = ", i))
      ## the next two lines really should be put outside the loop
      sigma.i <- se.psi.hat[i] # sigma_i in the formulas
      psi.hat.i <- psi.hat[i]
      ## The conditional distribution of psi_i is a mixture with several components
      ## Refer to equation 3.5
      ## Term1 comprises a Normal component restricted to (-infinity, psi)
      ## Term2 comprises a Normal component restricted to (psi, infinity)
      ## Term3 comprises a mixture of singletons (psi_j < psi)
      ## Term4 comprises a mixture of singletons (psi_j > psi)
      ## We need to first normalize everything so that we know which component of the
      ## mixture to pick.
      ## So Begin Normaliziing
      values <- psi.vec.current[-i]
      ## Calculate term3 of eqn 3.5
      term3.indices <- which (values < psi.current)
      m.minus <- length(term3.indices)
      psi.hat.term3 <- numeric(0) # don't need (never used)
      term3.summands <- numeric(0) # don't need
      term3.sum <- 0 # don't need
      if (m.minus > 0) {
        psi.vec.term3 <- values[term3.indices]
        term3.summands <- dnorm(psi.vec.term3, mean = psi.hat.i, sd = sigma.i) /
          (0.5 * M + m.minus)
        term3.sum <- sum(term3.summands)
      }
      ## print(paste("M minus is ", m.minus, " term3 sum", term3.sum))
      ##
      ## Calculate term4 of eqn 3.5
      term4.indices <- which (values > psi.current)
      m.plus <- length(term4.indices)
      psi.vec.term4 <- numeric(0) # don't need
      term4.summands <- numeric(0) # don't need
      term4.sum <- 0 # don't need
      if (m.plus > 0) {
        psi.vec.term4 <- values[term4.indices]
        term4.summands <- dnorm(psi.vec.term4, mean = psi.hat.i, sd = sigma.i) /
          (0.5 * M + m.plus)
        term4.sum <- sum(term4.summands)
      }
      ## print(paste("M plus is ", m.plus, " term4 sum ", term4.sum))
      ##
      ## Remember that despite the appearance, the above terms are mixtures
      ## of distributions we you need to keep them as vectors because later on,
      ## we'll have to decide which component of the mixture to pick, based
      ## on a uniform.
      ##
      ## print(paste("psi.current", psi.current, "tau.current", tau.current,
      ##              "sigma.i", sigma.i, "psi.hat.i", psi.hat.i))
      capital.A <- (psi.current * sigma.i^2 + psi.hat.i * tau.current^2) /
        (sigma.i^2 + tau.current^2)
      capital.B <- sqrt((sigma.i^2 * tau.current^2) / (sigma.i^2 + tau.current^2))
      common.factor <- dnorm(psi.current, mean = psi.hat.i,
                             sd = sqrt(sigma.i^2 + tau.current^2))
      ## print(paste("Capital A", capital.A, "Capital B", capital.B, "Common factor", common.factor))
      c.minus <- M * common.factor / (0.5 * M + m.minus)
      c.plus <- M * common.factor / (0.5 * M + m.plus)
      normalizing.constant <- (c.minus - c.plus) *
        pnorm(psi.current, mean = capital.A, sd = capital.B) +
          c.plus + term3.sum + term4.sum
      subdistribution.factor <- pnorm(psi.current, mean = capital.A, sd = capital.B)
      c.minus <- c.minus * subdistribution.factor
      c.plus <- c.plus * (1 - subdistribution.factor)
      ##    As a check, (c.minus + c.plus + term3.sum + term4.sum) / normalizing.constant = 1
      ##
      ##       print(paste("Check that ",
      ##                   c.minus, " + ", c.plus, " + ", term3.sum, " + ", term4.sum,
      ##                   " EQUALS ",
      ##                   normalizing.constant,
      ##                   (c.minus + c.plus + term3.sum + term4.sum) / normalizing.constant == 1))
      ##
      ## End of Normalization

      ## Now pick a mixture component randomly according to the mixing probabilities
      ## print(paste("C. minus", c.minus))
      ## print(paste("C. plus", c.plus))
      ## print(paste("term3 summands", term3.summands))
      ## print(paste("term4 summands", term4.summands))
      probability.vector <- c(c.minus, c.plus, term3.summands, term4.summands)
      ##         print(paste("Probability Vector", probability.vector))
      ##         print(probability.vector/normalizing.constant)
      ##
      ## In the sample call below, I start from -1 because I can use -1 and 0 to denote
      ## the Normal(-infinity, psi) and Normal(psi, infinity) components respy.
      mixture.component.index <- sample(-1:(m.minus+m.plus),
                                        1,
                                        replace=FALSE,
                                        prob = probability.vector)
      ## print(paste("Mixture component Index", mixture.component.index))
      a <- pnorm(psi.current, mean=capital.A, sd=capital.B)
      if (mixture.component.index <= -1) { # The Normal(-infinity, psi)
        psi.vec.current[i] <- capital.B * qnorm(a * runif(n=1)) + capital.A
      } else if (mixture.component.index <= 0) { # The Normal(psi, infinity)
        psi.vec.current[i] <- capital.B * qnorm(a + (1 - a) * runif(n=1)) + capital.A
      } else if (mixture.component.index <= m.minus) { # The term3 atoms
        psi.vec.current[i] <- values[term3.indices[mixture.component.index]]
      } else { # the term4 atoms
        psi.vec.current[i] <- values[term4.indices[mixture.component.index - m.minus]]
      }
    }
    psi.vec.current.unique <- unique(psi.vec.current)
    ## print(paste("Unique Psi", psi.vec.current.unique))
    m.star <- length(psi.vec.current.unique)
    ## print(paste("M star is ", m.star, "IC = ", ic, "I", i))
    aa.prime <- aa + 0.5 * m.star
    psi.vec.current.unique.bar.star <- mean(psi.vec.current.unique)
    bb.prime <- bb + 0.5 * (m.star * (psi.vec.current.unique.bar.star - cc)^2) /
      (1.0 + m.star * dd)
    if (length(psi.vec.current.unique) > 1) {
      bb.prime <- bb.prime + 0.5 * (m.star - 1) * var(psi.vec.current.unique)
    } # I don't see the point of this construct
    cc.prime <- (cc + m.star * dd * psi.vec.current.unique.bar.star) /
      (m.star * dd + 1.0)
    dd.prime <- dd / (1.0 + m.star * dd)
    dprime <- c(aa.prime,bb.prime,cc.prime,dd.prime)
    ## print(paste("aa.prime", aa.prime, "bb.prime", bb.prime, "psi.current", psi.current, "cc.prime", cc.prime, "dd.prime", dd.prime))
    psi.current <- rt.star.k(psi.vec.current.unique, m.star, M, dprime)
    tau.current <- 1/sqrt(rgamma(1, aa.prime + 0.5) /
                          (bb.prime + 0.5 * (psi.current - cc.prime)^2 / dd.prime))
    cycle.number <- cycle.number + 1 #20110701 deleted from output XXX
    output.matrix[ic+1,] <- c(psi.vec.current,
                              psi.current, tau.current)
  }
  z <- list(call = cl, ncycles=ncycles, M=M, prior = d,
            chain = output.matrix,start.user=start.user,start=start)
  class(z) <- c("dir.cond")
  z
}
