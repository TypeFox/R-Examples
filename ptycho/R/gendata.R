createData <- function(X, y, omega=NULL, beta=NULL) {
  if (!is.matrix(X)) {
    if (!is.list(X) || length(X) != 2) {
      stop("createData: X must be matrix or list of length 2")
    }
    X <- do.call(X[[1]], X[[2]])
  }
  p <- ncol(X)
  if (is.numeric(y)) {
    if (is.vector(y)) dim(y) <- c(nrow(X), 1)
    if (nrow(y) != nrow(X)) {
      stop("createData: X and y must have same number of rows")
    }
    q <- ncol(y)
    noise.sd <- NULL
    eta2 <- computeEta2(X, y)
    reps <- list(list(omega=NULL, omega.grp=NULL, indic.grp=NULL,
                      indic.var=NULL, tau=NULL, beta=NULL, y=y, eta2=eta2))
  } else if (is.list(y)) {
    if (any(is.na(match(c("nreps","q","sd"), names(y))))) {
      stop(gettextf("createData: list y must have tags %s, %s, %s",
                    dQuote("nreps"), dQuote("q"), dQuote("sd")))
    }
    q <- y$q
    noise.sd <- y$sd
    reps <- simulateReplicates(X, omega, beta, y)
  } else {
    stop("createData: y must be vector, matrix, or list")
  }
  list(X=X, q=q, noise.sd=noise.sd, omega=omega, beta=beta, replicates=reps)
}

createDataBayesModel <- function(mode=c("exchange","pleiotropy","gene"),
                                 n, p, q, nreps, tau.min, tau.max, G=NULL) {
  mode <- match.arg(mode)
  if (mode == "exchange") {
    omegafn <- "createOmegaBeta"
    omegaargs <- list(shape1=12, shape2=48, byrow=TRUE)
  } else if (mode == "pleiotropy") {
    omegafn <- "createOmegaCrossTraits"
    omegaargs <- list(indic.grp.shape1=16, indic.grp.shape2=55,
                      shape1=48, shape2=12)
  } else if (mode=="gene") {
    groups <- createGroupsSim(G, p)
    omegafn <- "createOmegaCrossVars"
    omegaargs <- list(indic.grp.shape1=16, indic.grp.shape2=55,
                      shape1=48, shape2=12, groups=groups)
  }
  data <- createData(X=list("createOrthogonalX", list(n=n,p=p)),
                     y=list(nreps=nreps, q=q, sd=1),
                     omega=list(omegafn, omegaargs),
                     beta=list("createBetaNormal",
                               list(tau.min=tau.min, tau.max=tau.max)))
}

drawIndicator <- function(nrow, ncol, prob) {
  n <- nrow*ncol
  repeat {
    x <- rbinom(n=n, size=1, prob=prob)
    dim(x) <- c(nrow, ncol)
    if (all(colSums(x) > 0)) return(x)
    #cat("REDRAWING\n")
  }
}

simulateReplicates <- function(X, omega, beta, y) {
  z <- rlply(y$nreps, simulateReplicate(X, omega, beta, y))
}

simulateReplicate <- function(X, omega, beta, y) {
  p <- ncol(X)
  q <- y$q
  noise.sd <- y$sd
  # The column names of the responses should also be put on beta and such.
  ynames <- createYNames(q)
  # Generate omega
  z <- createOmega(p, q, omega)
  # The column names of X may or may not be bogus.
  rownames(z$omega) <- colnames(X)
  colnames(z$omega) <- ynames
  # Generate indic.var
  indic.var <- drawIndicator(nrow=p, ncol=q, prob=z$omega)
  rownames(indic.var) <- colnames(X)
  colnames(indic.var) <- ynames
  z <- append(z, list(indic.var=indic.var))
  # Generate beta
  z <- append(z, createBeta(z$indic.var, noise.sd, beta))
  # Generate y
  z <- append(z, createY(X, z$beta, noise.sd, ynames))
  z
}

createOmega <- function(p, q, omega) {
  if (is.list(omega)) {
    if (length(omega) != 2) {
      stop("createData: list omega must have length 2")
    }
    z <- do.call(omega[[1]], append(list(p=p,q=q), omega[[2]]))
  } else {
    z <- augmentOmega(omega, p, q)
  }
  z
}

augmentOmega <- function(omega, p, q) {
  if (is.matrix(omega)) {
    if (nrow(omega) != p || ncol(omega) != q) {
      stop("createData: omega must be p-by-q")
    }
  } else {
    if (length(omega) != 1 && length(omega) != p) {
      warning("createData: omega is a weird size")
    }
    omega <- matrix(omega, nrow=p, ncol=q)
  }
  list(omega=omega, omega.grp=NULL, indic.grp=NULL)
}

# Beta in the name of this function refers to the Beta distribution (not to the
# covariate effect sizes).  The returned values are generated directly from the
# beta distribution specified by the function arguments.
createOmegaBeta <- function(p, q, shape1, shape2, byrow) {
  if (byrow) {
    ngen <- q; nmult <- p
  } else {
    ngen <- p; nmult <- q
  }
  z <- rbeta(ngen, shape1, shape2)
  z <- matrix(z, nrow=p, ncol=q, byrow=byrow)
  list(omega=z, omega.grp=NULL, indic.grp=NULL)
}

# Generate omega as in model with multiple responses and doGrpIndicator
createOmegaCrossTraits <- function(p, q, indic.grp.shape1, indic.grp.shape2,
                                   shape1, shape2) {
  omega.grp <- rbeta(1, indic.grp.shape1, indic.grp.shape2)
  omega.grp <- matrix(omega.grp, nrow=p, ncol=1)
  rownames(omega.grp) <- createXNames(p)
  colnames(omega.grp) <- "group"
  indic.grp <- drawIndicator(nrow=p, ncol=1, prob=omega.grp)
  rownames(indic.grp) <- rownames(omega.grp)
  colnames(indic.grp) <- colnames(omega.grp)
  omega <- rbeta(p, shape1, shape2) * indic.grp
  omega <- matrix(omega, nrow=p, ncol=q)
  list(omega=omega, omega.grp=omega.grp, indic.grp=indic.grp)
}

# Generate omega as in model that combines across covariates and
# doGrpIndicator=T.
createOmegaCrossVars <- function(p, q, indic.grp.shape1, indic.grp.shape2,
                                 shape1, shape2, groups) {
  G <- length(groups$group2var)
  omega.grp <- rbeta(q, indic.grp.shape1, indic.grp.shape2)
  omega.grp <- matrix(omega.grp, nrow=G, ncol=q, byrow=TRUE)
  rownames(omega.grp) <- names(groups$group2var)
  colnames(omega.grp) <- createYNames(q)
  indic.grp <- drawIndicator(nrow=G, ncol=q, prob=omega.grp)
  rownames(indic.grp) <- rownames(omega.grp)
  colnames(indic.grp) <- colnames(omega.grp)
  omega <- rbeta(p*q, shape1, shape2)
  dim(omega) <- c(p, q)
  omega <- omega * indic.grp[groups$var2group,]
  list(omega=omega, omega.grp=omega.grp, indic.grp=indic.grp)
}

# Generate omega when using actual genotype data.
createOmegaActualX <- function(p, q, n.common, probs.common, groups) {
  if (length(n.common) != length(probs.common)) {
    stop("createData: n.common and probs.common must have same length")
  }
  # Initialize indic.grp and omega
  indic.grp <- matrix(0, nrow=max(groups$var2group), ncol=q)
  gnames <- names(groups$group2var)
  if (length(gnames) < nrow(indic.grp)) {
    # Some groups have been omitted because they should not be selected for
    # inclusion in a true model.
    varmsg <- which(groups$var2group > length(gnames))
    gnames <- c(gnames, names(groups$var2group)[varmsg])
  }
  rownames(indic.grp) <- gnames
  colnames(indic.grp) <- createYNames(q)
  omega <- matrix(0, nrow=p, ncol=q)
  # For each response, select a gene with at least five rare variants.
  grp.size <- laply(groups$group2var, length)
  grpsel <- sample(which(grp.size >= 5), size=q, replace=TRUE)
  # Select rare variants
  for (nn in seq_len(q)) {
    indic.grp[grpsel[nn],nn] <- 1
    nvar <- sample(3:4, size=1)
    vars <- sample(groups$group2var[[grpsel[nn]]], size=nvar)
    omega[vars,nn] <- 1
  }
  # Select common variants
  # Must select all at the same time to ensure no duplicates
  grpsel <- sample(which(grp.size == 1), size=sum(n.common))
  ncum <- c(0, cumsum(n.common))
  for (nn in seq_along(n.common)) {
    if (ncum[nn] < ncum[nn+1]) {
      grptmp <- grpsel[(ncum[nn]+1):ncum[nn+1]]
      indic.grp[grptmp,] <- 1
      for (ngrp in grptmp) omega[groups$group2var[[ngrp]],] <- probs.common[nn]
    }
  }
  list(omega=omega, omega.grp=NULL, indic.grp=indic.grp)
}

# Generate omega for test using correlated variables
make.createOmegaCorTest <- function(nreps) {
  count <- 0
  f <- function(p, q, vars) {
    if (length(vars) != 2) {
      stop("createOmegaCorTest: must specify two variants")
    }
    count <<- count + 1
    omega <- matrix(0, nrow=p, ncol=q)
    if (count <= nreps/2) {
      omega[vars[1],1] <- 1
      omega[vars[2],2] <- 1
    } else if (count <= 3*nreps/4) {
      omega[vars[1],] <- 1
    } else {
      omega[vars[2],] <- 1
    }
    list(omega=omega, omega.grp=NULL, indic.grp=NULL)
  }
  return(f)
}

createBeta <- function(indic.var, noise.sd, beta) {
  if (!is.list(beta) || length(beta) != 2) {
    stop("createData: beta must be list of length 2")
  }
  z <- do.call(beta[[1]],
               append(list(indic.var=indic.var, noise.sd=noise.sd), beta[[2]]))
  rownames(z$beta) <- rownames(indic.var)
  colnames(z$beta) <- colnames(indic.var)
  z
}

# For indic.var=1, \abs{beta} is uniform.
createBetaUnif <- function(indic.var, noise.sd, min, max) {
  p <- nrow(indic.var); q <- ncol(indic.var)
  beta <- sample(c(-1,1), size=p*q, replace=TRUE)
  beta <- beta * runif(p*q, min=min, max=max)
  beta <- beta * indic.var
  list(tau=NULL, beta=beta)
}

# For indic.var=1, beta is normal.
createBetaNormal <- function(indic.var, noise.sd, tau.min, tau.max) {
  p <- nrow(indic.var); q <- ncol(indic.var)
  tau <- runif(1, min=tau.min, max=tau.max)
  beta <- rnorm(p*q, sd=noise.sd*tau)
  beta <- beta * indic.var
  list(tau=tau, beta=beta)
}

createY <- function(X, beta, noise.sd, ynames) {
  y <- X %*% beta
  ydims <- dim(y)
  noise <- rnorm(prod(ydims), sd=noise.sd)
  dim(noise) <- ydims
  y <- y + noise
  colnames(y) <- ynames
  if (max(abs(colSums(X))) < 1e-6) y <- scale(y, scale=FALSE)
  eta2 <- computeEta2(X, y)
  list(y=y, eta2=eta2)
}

computeEta2 <- function(X, y) {
  Xty <- t(X) %*% y
  yty <- aaply(y, .margins=2, function(z) { sum(z^2) })
  eta2 <- scale(Xty^2, center=FALSE, scale=nrow(y)*yty)
  rownames(eta2) <- colnames(X)
  colnames(eta2) <- colnames(y)
  eta2
}

createPubData <- function(mode=c("tinysim","ptychoIn",
                                 "exchange","pleiotropy","gene",
                                 "actualGeno","actualPheno","corTest",
                                 "fixedOmega","uniformEffects"),
                          X=NULL, y=NULL, var.detail=NULL, variants=NULL) {
  mode <- match.arg(mode)
  if (mode != "corTest") set.seed(1234)
  if (mode == "actualPheno") {
    # Actual data
    data <- createData(X, y)
  } else if (mode == "actualGeno") {
    # Using actual genotypes, simulate phenotypes as in paper
    groups <- createGroupsActualX(X, var.detail,
                                  MAF.threshold=0.01, is.sim=TRUE)
    data <- createData(X, y=list(nreps=100, q=3, sd=1),
                       omega=list("createOmegaActualX",
                                  list(n.common=c(10,40),
                                       probs.common=c(0.9,0.1), groups=groups)),
                       beta=list("createBetaNormal",
                                 list(tau.min=0.045, tau.max=0.063)))
  } else if (mode == "corTest") {
    nreps <- 40
    createOmegaCorTest <- make.createOmegaCorTest(nreps=nreps)
    data <- createData(X, y=list(nreps=nreps, q=2, sd=1),
                       omega=list(createOmegaCorTest,
                                  list(vars=which(colnames(X) %in% variants))),
                       beta=list("createBetaNormal",
                                 list(tau.min=0.045, tau.max=0.063)))
  } else if (mode == "fixedOmega") {
    # orthogonal X, fixed omega
    data <- createData(X=list("createOrthogonalX", list(n=5000,p=50)),
                       y=list(nreps=100, q=5, sd=1),
                       omega=c(rep(0.9,10), rep(0.1,40)),
                       beta=list("createBetaNormal",
                                 list(tau.min=0.045, tau.max=0.063)))
  } else if (mode == "uniformEffects") {
    # orthogonal X, uniform effect sizes
    data <- createData(X=list("createOrthogonalX", list(n=5000,p=50)),
                       y=list(nreps=100, q=5, sd=1),
                       omega=c(rep(0.9,10), rep(0.1,40)),
                       beta=list("createBetaUnif", list(min=0.01,max=0.1)))
  } else {
    if (mode == "tinysim") {
      mode <- "pleiotropy"
      n <- 100; p <- 10; q <- 5; nreps <- 10
    } else if (mode == "ptychoIn") {
      set.seed(4)
      mode <- "gene"
      n <- 3000; p <- 10; q <- 1; nreps <- 1
    } else {
      n <- 5000; p <- 50; q <- 5; nreps <- 100
    }
    G <- if (mode=="gene") p/5 else NULL
    data <- createDataBayesModel(mode, n, p, q, nreps, 
                                 tau.min=0.045, tau.max=0.063, G=G)
  }
}
