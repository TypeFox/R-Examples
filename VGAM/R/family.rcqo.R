# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.






rcqo <- function(n, p, S,
                 Rank = 1,
                 family = c("poisson", "negbinomial", "binomial-poisson",
                            "Binomial-negbinomial", "ordinal-poisson",
                            "Ordinal-negbinomial", "gamma2"),
                 eq.maximums = FALSE,
                 eq.tolerances = TRUE,
                 es.optimums = FALSE,
                 lo.abundance = if (eq.maximums) hi.abundance else 10,
                 hi.abundance = 100,
                 sd.latvar = head(1.5/2^(0:3), Rank),
                 sd.optimums = ifelse(es.optimums, 1.5/Rank, 1) *
                            ifelse(scale.latvar, sd.latvar, 1),
                 sd.tolerances = 0.25,
                 Kvector = 1,
                 Shape = 1,
                 sqrt.arg = FALSE,
                 log.arg = FALSE,
                 rhox = 0.5,
                 breaks = 4,  # ignored unless family = "ordinal"
                 seed = NULL,
                 optimums1.arg = NULL,
                 Crow1positive = TRUE,
                 xmat = NULL,  # Can be input
                 scale.latvar = TRUE) {
  family <- match.arg(family,
                      c("poisson", "negbinomial", "binomial-poisson",
                        "Binomial-negbinomial", "ordinal-poisson",
                        "Ordinal-negbinomial", "gamma2"))[1]

  if (!is.Numeric(n, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1))
    stop("bad input for argument 'n'")
  if (!is.Numeric(p, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1) ||
      p < 1 + Rank)
    stop("bad input for argument 'p'")
  if (!is.Numeric(S, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1))
    stop("bad input for argument 'S'")
  if (!is.Numeric(Rank, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1) ||
      Rank > 4)
    stop("bad input for argument 'Rank'")
  if (!is.Numeric(Kvector, positive = TRUE))
    stop("bad input for argument 'Kvector'")
  if (!is.Numeric(rhox) || abs(rhox) >= 1)
    stop("bad input for argument 'rhox'")
  if (length(seed) &&
      !is.Numeric(seed, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'seed'")
  if (!is.logical(eq.tolerances) ||
      length(eq.tolerances) > 1)
    stop("bad input for argument 'eq.tolerances)'")
  if (!is.logical(sqrt.arg) || length(sqrt.arg) > 1)
    stop("bad input for argument 'sqrt.arg)'")
  if (family != "negbinomial" && sqrt.arg)
    warning("argument 'sqrt.arg' is used only with family='negbinomial'")
  if (!eq.tolerances && !is.Numeric(sd.tolerances, positive = TRUE))
    stop("bad input for argument 'sd.tolerances'")
  if (!is.Numeric(lo.abundance, positive = TRUE))
    stop("bad input for argument 'lo.abundance'")
  if (!is.Numeric(sd.latvar, positive = TRUE))
    stop("bad input for argument 'sd.latvar'")
  if (!is.Numeric(sd.optimums, positive = TRUE))
    stop("bad input for argument 'sd.optimums'")
  if (eq.maximums && lo.abundance != hi.abundance)
    stop("arguments 'lo.abundance' and 'hi.abundance' must ",
         "be equal when 'eq.tolerances = TRUE'")
  if (any(lo.abundance > hi.abundance))
    stop("lo.abundance > hi.abundance is not allowed")
  if (!is.logical(Crow1positive)) {
    stop("bad input for argument 'Crow1positive)'")
  } else {
    Crow1positive <- rep(Crow1positive, len = Rank)
  }
  Shape <- rep(Shape, len = S)
  sd.latvar <- rep(sd.latvar, len = Rank)
  sd.optimums <- rep(sd.optimums, len = Rank)
  sd.tolerances <- rep(sd.tolerances, len = Rank)
  AA <- sd.optimums / 3^0.5
  if (Rank > 1 && any(diff(sd.latvar) > 0))
   stop("argument 'sd.latvar)' must be a vector with decreasing values")

  if (FALSE)
  change.seed.expression <- expression({
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      runif(1)  # initialize the RNG if necessary
    }
    if (is.null(seed)) {
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
  })


  change.seed.expression <- expression({
      if (length(seed)) set.seed(seed)
  })
  eval(change.seed.expression)

  V <- matrix(rhox, p-1, p-1)
  diag(V) <- 1
  L <- chol(V)
  if (length(xmat)) {
    xnames <- colnames(xmat)
  } else {
    eval(change.seed.expression)
    xmat <- matrix(rnorm(n*(p-1)), n, p-1) %*% L
    xmat <- scale(xmat, center = TRUE)
    xnames <- paste("x", 2:p, sep = "")
    dimnames(xmat) <- list(as.character(1:n), xnames)
  }
  eval(change.seed.expression)
  Ccoefs <- matrix(rnorm((p-1)*Rank), p-1, Rank)
  latvarmat <- cbind(xmat %*% Ccoefs)
  if (Rank > 1) {
    Rmat <- chol(var(latvarmat))
    iRmat <- solve(Rmat)
    latvarmat <- latvarmat %*% iRmat  # var(latvarmat) == diag(Rank)
    Ccoefs <- Ccoefs %*% iRmat
  }
  for (r in 1:Rank)
    if (( Crow1positive[r] && Ccoefs[1, r] < 0) ||
        (!Crow1positive[r] && Ccoefs[1, r] > 0)) {
      Ccoefs[ , r] <- -Ccoefs[ , r]
      latvarmat[ , r] <- -latvarmat[ , r]
    }

  if (scale.latvar) {
    for (r in 1:Rank) {
      sd.latvarr <- sd(latvarmat[, r])
      latvarmat[, r] <- latvarmat[, r] * sd.latvar[r] / sd.latvarr
      Ccoefs[, r]  <- Ccoefs[, r] * sd.latvar[r] / sd.latvarr
    }
  } else {
    sd.latvarr <- NULL
    for (r in 1:Rank) {
      sd.latvarr <- c(sd.latvarr, sd(latvarmat[, r]))
    }
  }
  if (es.optimums) {
    if (!is.Numeric(S^(1/Rank), integer.valued = TRUE) ||
        S^(1/Rank) < 2)
      stop("S^(1/Rank) must be an integer greater or equal to 2")
    if (Rank == 1) {
      optimums <- matrix(NA_real_, S, Rank)
      for (r in 1:Rank) {
        optimums[, r] <- seq(-AA, AA, len = S^(1/Rank))
      }
    } else if (Rank == 2) {
      optimums <-
        expand.grid(latvar1 = seq(-AA[1], AA[1], len = S^(1/Rank)),
                    latvar2 = seq(-AA[2], AA[2], len = S^(1/Rank)))
    } else if (Rank == 3) {
      optimums <-
        expand.grid(latvar1 = seq(-AA[1], AA[1], len = S^(1/Rank)),
                    latvar2 = seq(-AA[2], AA[2], len = S^(1/Rank)),
                    latvar3 = seq(-AA[3], AA[3], len = S^(1/Rank)))
    } else {
      optimums <-
        expand.grid(latvar1 = seq(-AA[1], AA[1], len = S^(1/Rank)),
                    latvar2 = seq(-AA[2], AA[2], len = S^(1/Rank)),
                    latvar3 = seq(-AA[3], AA[3], len = S^(1/Rank)),
                    latvar4 = seq(-AA[4], AA[4], len = S^(1/Rank)))
    }
    if (Rank > 1)
    optimums <- matrix(unlist(optimums), S, Rank)  # Make sure its a matrix
  } else {
    optimums <- matrix(1, S, Rank)
    eval(change.seed.expression)
    for (r in 1:Rank) {
      optimums[, r] <- rnorm(n = S, sd = sd.optimums[r])
    }
  }
  for (r in 1:Rank)
    optimums[, r] <- optimums[, r] * sd.optimums[r] / sd(optimums[, r])


  if (length(optimums1.arg) && Rank == 1)
  for (r in 1:Rank)
    optimums[, r] <- optimums1.arg



  ynames <- paste("y", 1:S, sep = "")
  Kvector <- rep(Kvector, len = S)
  names(Kvector) <- ynames
  latvarnames <- if (Rank == 1) "latvar" else
                 paste("latvar", 1:Rank, sep = "")
  Tols <- if (eq.tolerances) {
    matrix(1, S, Rank)
  } else {
    eval(change.seed.expression)
    temp <- matrix(1, S, Rank)
    if (S > 1)
      for (r in 1:Rank) {
        temp[-1, r] <- rnorm(S-1, mean = 1, sd = sd.tolerances[r])
        if (any(temp[, r] <= 0))
          stop("negative tolerances!")
        temp[, r] <- temp[, r]^2  # Tolerance matrix  = var-cov matrix)
      }
    temp
  }

  dimnames(Tols)   <- list(ynames, latvarnames)
  dimnames(Ccoefs) <- list(xnames, latvarnames)
  dimnames(optimums) <- list(ynames, latvarnames)
  loeta <- log(lo.abundance)   # May be a vector
  hieta <- log(hi.abundance)
  eval(change.seed.expression)
  log.maximums <- runif(S, min = loeta, max = hieta)
  names(log.maximums) <- ynames
  etamat <- matrix(log.maximums, n, S, byrow = TRUE)
  for (jay in 1:S) {
    optmat <- matrix(optimums[jay, ], nrow = n, ncol = Rank, byrow = TRUE)
    tolmat <- matrix(    Tols[jay, ], nrow = n, ncol = Rank, byrow = TRUE)
    temp <- cbind((latvarmat - optmat) / tolmat)
    for (r in 1:Rank)
      etamat[, jay] <- etamat[, jay] -
                       0.5 * (latvarmat[, r] - optmat[jay, r]) * temp[, r]
  }

  rootdist <- switch(family,
    "poisson" = 1, "binomial-poisson" = 1, "ordinal-poisson" = 1,
    "negbinomial" = 2, "Binomial-negbinomial" = 2,
    "Ordinal-negbinomial" = 2,
    "gamma2" = 3)
  eval(change.seed.expression)
  if (rootdist == 1) {
    ymat <- matrix(rpois(n * S, lambda = exp(etamat)), n, S)
  } else if (rootdist == 2) {
    mKvector <- matrix(Kvector, n, S, byrow = TRUE)
    ymat <- matrix(rnbinom(n = n * S, mu = exp(etamat), size = mKvector),
                   n, S)
    if (sqrt.arg)
      ymat <- ymat^0.5
  } else if (rootdist == 3) {
    Shape <- matrix(Shape, n, S, byrow = TRUE)
    ymat <- matrix(rgamma(n * S, shape = Shape,
                                 scale = exp(etamat) / Shape),
                   n, S)
    if (log.arg) ymat <- log(ymat)
  } else {
    stop("argument 'rootdist' unmatched")
  }

  tmp1 <- NULL
  if (any(family == c("ordinal-poisson", "Ordinal-negbinomial"))) {
    tmp1 <- cut(c(ymat), breaks = breaks, labels = NULL)
    ymat <- cut(c(ymat), breaks = breaks, labels = FALSE)
    dim(ymat) <- c(n,S)
    }
    if (any(family == c("binomial-poisson", "Binomial-negbinomial")))
      ymat <- 0 + (ymat > 0)

    myform <- as.formula(paste(paste("cbind(",
             paste(paste("y", 1:S, sep = ""), collapse = ", "),
             ") ~ ", sep = ""),
             paste(paste("x", 2:p, sep = ""), collapse = "+"), sep = ""))

  dimnames(ymat) <- list(as.character(1:n), ynames)
  ans <- data.frame(xmat, ymat)
  attr(ans, "concoefficients") <- Ccoefs
  attr(ans, "Crow1positive") <- Crow1positive
  attr(ans, "family") <- family
  attr(ans, "formula") <- myform  # Useful for running cqo() on the data
  attr(ans, "Rank") <- Rank
  attr(ans, "family") <- family
  attr(ans, "Kvector") <- Kvector
  attr(ans, "log.maximums") <- log.maximums
  attr(ans, "lo.abundance") <- lo.abundance
  attr(ans, "hi.abundance") <- hi.abundance
  attr(ans, "optimums") <- optimums
  attr(ans, "log.arg") <- log.arg
  attr(ans, "latvar") <- latvarmat
  attr(ans, "eta") <- etamat
  attr(ans, "eq.tolerances") <- eq.tolerances
  attr(ans, "eq.maximums") <- eq.maximums ||
                              all(lo.abundance == hi.abundance)
  attr(ans, "es.optimums") <- es.optimums
  attr(ans, "seed") <- seed  # RNGstate
  attr(ans, "sd.tolerances") <- sd.tolerances
  attr(ans, "sd.latvar") <- if (scale.latvar) sd.latvar else sd.latvarr
  attr(ans, "sd.optimums") <- sd.optimums
  attr(ans, "Shape") <- Shape
  attr(ans, "sqrt") <- sqrt.arg
  attr(ans, "tolerances") <- Tols^0.5  # Like a standard deviation
  attr(ans, "breaks") <- if (length(tmp1)) attributes(tmp1) else breaks
  ans
}




 if (FALSE)
dcqo <-
  function(x, p, S,
           family = c("poisson", "binomial", "negbinomial", "ordinal"),
           Rank = 1,
           eq.tolerances = TRUE,
           eq.maximums = FALSE,
           EquallySpacedOptima = FALSE,
           lo.abundance = if (eq.maximums) 100 else 10,
           hi.abundance = 100,
           sd.tolerances = 1,
           sd.optimums = 1,
           nlevels = 4,  # ignored unless family = "ordinal"
           seed = NULL) {
 warning("20060612; needs a lot of work based on rcqo()")


  if (mode(family) != "character" && mode(family) != "name")
    family <- as.character(substitute(family))
  family <- match.arg(family, c("poisson", "binomial",
                                 "negbinomial", "ordinal"))[1]


  if (!is.Numeric(p, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1) ||
      p < 2)
    stop("bad input for argument 'p'")
  if (!is.Numeric(S, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1))
    stop("bad input for argument 'S'")
  if (!is.Numeric(Rank, integer.valued = TRUE,
                  positive = TRUE, length.arg = 1))
    stop("bad input for argument 'Rank'")
  if (length(seed) &&
      !is.Numeric(seed, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'seed'")
  if (!is.logical(eq.tolerances) || length(eq.tolerances)>1)
    stop("bad input for argument 'eq.tolerances)'")
  if (eq.maximums && lo.abundance != hi.abundance)
    stop("'lo.abundance' and 'hi.abundance' must ",
         "be equal when 'eq.tolerances = TRUE'")
  if (length(seed)) set.seed(seed)

  xmat <- matrix(rnorm(n*(p-1)), n, p-1,
                 dimnames = list(as.character(1:n),
                                 paste("x", 2:p, sep = "")))
  Ccoefs <- matrix(rnorm((p-1)*Rank), p-1, Rank)
  latvarmat <- xmat %*% Ccoefs
  optimums <- matrix(rnorm(Rank*S, sd = sd.optimums), S, Rank)
  Tols <- if (eq.tolerances) matrix(1, S, Rank) else
         matrix(rnorm(Rank*S, mean = 1, sd = 1), S, Rank)
  loeta <- log(lo.abundance)
  hieta <- log(hi.abundance)
  log.maximums <- runif(S, min = loeta, max = hieta)

  etamat <- matrix(log.maximums, n, S, byrow = TRUE)
  for (jay in 1:S) {
    optmat <- matrix(optimums[jay, ], n, Rank, byrow = TRUE)
    tolmat <- matrix(  Tols[jay, ], n, Rank, byrow = TRUE)
    temp <- cbind((latvarmat - optmat) * tolmat)
    for (r in 1:Rank)
      etamat[, jay] <- etamat[, jay] - 0.5 * temp[, r] *
                       (latvarmat[, r] - optmat[jay, r])
  }

  ymat <- if (family == "negbinomial") {



  } else {
     matrix(rpois(n * S, lambda = exp(etamat)), n, S)
  }
  if (family == "binomial")
    ymat <- 0 + (ymat > 0)

  dimnames(ymat) <- list(as.character(1:n),
                         paste("y", 1:S, sep = ""))
  ans <- data.frame(xmat, ymat)
  attr(ans, "concoefficients") <- Ccoefs
  attr(ans, "family") <- family
  ans
}





getInitVals <- function(gvals, llfun, ...) {
  LLFUN <- match.fun(llfun)
  ff <- function(myx, ...) LLFUN(myx, ...)
  objFun <- gvals
  for (ii in 1:length(gvals))
    objFun[ii] <- ff(myx = gvals[ii], ...) 
  try.this <- gvals[objFun == max(objFun)]  # Usually scalar, maybe vector
  try.this
}









campp <- function(q, size, prob, mu) {
  if (!missing(mu)) {
    if (!missing(prob))
      stop("arguments 'prob' and 'mu' both specified")
    prob <- size/(size + mu)
  }
  K <- (1/3) * ((9*q+8) / (q+1) -
       ((9*size-1)/size) *
       (mu/(q+1))^(1/3)) / sqrt( (1/size) *
                   (mu/(q+1))^(2/3) + 1 / (q+1))  # Note the +, not -
  pnorm(K)
}










