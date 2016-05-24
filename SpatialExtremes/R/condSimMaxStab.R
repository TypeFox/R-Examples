condrmaxstab <- function(k = 1, coord, cond.coord, cond.data, cov.mod = "powexp",
                         ...,  do.sim = TRUE, thin = n.cond, burnin = 50, parts){

  if (!(cov.mod %in% c("brown", "whitmat", "powexp", "bessel", "cauchy",
                       "twhitmat", "tpowexp", "tbessel", "tcauchy")))
    stop("'cov.mod' must be one of 'brown', '(t)whitmat', '(t)powexp', '(t)bessel' or '(t)cauchy'")

  timings <- rep(NA, 3)
  if (is.null(dim(coord))){
    n.site <- length(coord)
    n.cond <- length(cond.coord)

    all.coord <- scale(c(cond.coord, coord), scale = FALSE)
    dim <- 1
  }

  else {
    n.site <- nrow(coord)
    n.cond <- nrow(cond.coord)

    all.coord <- scale(rbind(cond.coord, coord), scale = FALSE)
    dim <- ncol(coord)
  }

  dist <- as.matrix(dist(all.coord, diag = TRUE, upper = TRUE))

  n <- n.site + n.cond

  if (cov.mod == "brown")
    model <- "Brown-Resnick"

  else if (substr(cov.mod, 1, 1) == "t"){
      model <- "extremal-t"
      cov.mod <- substr(cov.mod, 2, 10)
  }

  else
      model <- "Schlather"

  ## Get the parameters
  if (model == "extremal-t")
      DoF <- list(...)$DoF

  range <- list(...)$range
  smooth <- list(...)$smooth

  if (is.null(range) || is.null(smooth))
    stop("You must specify values for the range and the smooth parameters")

  if (range <= 0)
    stop("'range' must be positive")

  if (model == "Brown-Resnick"){
    if ((smooth <= 0) || (smooth > 2))
      stop("smooth must belongs to (0, 2] for Brown-Resnick models")

    ## Compute some fixed quantities once for all
    one <- rep(1, n)
    one.sub <- rep(1, n.cond)

    if (dim == 1)
      sigma2 <- 2 * (abs(all.coord) / range)^smooth

    else {
      h <- sqrt(rowSums(all.coord^2))
      sigma2 <- 2 * (h / range)^smooth
    }

    sigma2.sub <- sigma2[1:n.cond]
    y <- log(cond.data)
    y.tilde <- c(y, rep(0, n.site))

    cov <- - (dist / range)^smooth
    for (i in 1:n)
      for (j in 1:n)
        cov[i,j] <- 0.5 * (sigma2[i] + sigma2[j]) + cov[i,j]

    cov.sub <- cov[1:n.cond, 1:n.cond]
    cov.chol <- chol(cov)
    cov.chol.sub <- chol(cov.sub)

    icov.chol1 <- backsolve(cov.chol,  one, transpose = TRUE)
    icov.chol1.sub <- backsolve(cov.chol.sub, one.sub, transpose = TRUE)
    imahal1 <- 1 / as.numeric(t(icov.chol1) %*% icov.chol1)
    imahal1.sub <- 1 / as.numeric(t(icov.chol1.sub) %*% icov.chol1.sub)
    ham <- diag(n) - icov.chol1 %*% t(icov.chol1) * imahal1
    ham.sub <- diag(n.cond) - icov.chol1.sub %*% t(icov.chol1.sub) * imahal1.sub
    icov.cholsigma2 <- backsolve(cov.chol, sigma2, transpose = TRUE)
    icov.cholsigma2.sub <- backsolve(cov.chol.sub, sigma2.sub, transpose = TRUE)
    mean1 <- -0.5 * icov.cholsigma2 + (0.5 * as.numeric(t(icov.chol1) %*% icov.cholsigma2) -
                                       1) * imahal1 * icov.chol1
    mean1.sub <- -0.5 * icov.cholsigma2.sub +
      (0.5 * as.numeric(t(icov.chol1.sub) %*% icov.cholsigma2.sub) - 1) *
        imahal1.sub * icov.chol1.sub
  }

  if (model %in% c("Schlather", "extremal-t")){
    y <- cond.data
    cov.fun <- covariance(nugget = 0, sill = 1, range = range, smooth = smooth,
                          cov.mod = cov.mod, plot = FALSE)

    cov <- cov.fun(dist)
    cov.sub <- cov[1:n.cond, 1:n.cond]
    cov.chol.sub <- chol(cov.sub)
  }

  if (missing(parts)){
    ## If n.cond > 7 we use a Gibbs sampler otherwise we compute all the
    ## weights and sample directly from their distribution

    if (n.cond <= 7){
      ## List all possible partitions
      start <- proc.time()
      n.part <- .C("bell", as.integer(n.cond), n.part = integer(1))$n.part
      all.part <- .C("listAllPartOfASet", as.integer(n.cond), as.integer(n.part),
                     all.part = integer(n.part * n.cond), all.size = integer(n.part))
      all.size <- all.part$all.size
      all.part <- all.part$all.part

      ## Compute the weight for each single partition

      if (model == "Brown-Resnick")
        weights <- .C("computeWeightsBR", as.integer(n.cond), as.double(y), as.integer(n.part),
                      as.integer(all.part), as.integer(all.size), as.double(cov.sub),
                      as.double(sigma2.sub), as.double(cov.chol.sub), as.double(ham.sub),
                      as.double(mean1.sub), weights = double(n.part))$weights

      if (model == "Schlather")
        weights <- .C("computeWeightsSC", as.integer(n.cond), as.double(y), as.integer(n.part),
                      as.integer(all.part), as.integer(all.size), as.double(cov.sub),
                      weights = double(n.part))$weights

      if (model == "extremal-t")
          weights <- .C("computeWeightsExtt", as.integer(n.cond), as.double(y), as.integer(n.part),
                      as.integer(all.part), as.integer(all.size), as.double(cov.sub), as.double(DoF),
                      weights = double(n.part))$weights

      all.part <- matrix(all.part, n.cond, n.part)
      idx.part <- sample(1:n.part, k, replace = TRUE, prob = weights)
      parts <- all.part[,idx.part]
      timings[1] <- (proc.time() - start)[3]
    }

    else {
      ## Get an appropriate starting partition by sampling max-stable
      ## processes and picking the most likely partition---quite naive
      ## though
      n.sim.start <- 250

      if (model == "Brown-Resnick")
        dummy <- .C("getStartingPartitionBR", as.integer(n.sim.start), as.integer(n.cond),
                    as.double(cond.coord), as.double(range), as.double(smooth),
                    start = integer(n.sim.start * n.cond))$start

      if (model == "Schlather")
        dummy <- .C("getStartingPartitionSC", as.integer(n.sim.start), as.integer(n.cond),
                    as.double(cov.chol.sub), start = integer(n.sim.start * n.cond))$start

      if (model == "extremal-t")
          dummy <- .C("getStartingPartitionExtt", as.integer(n.sim.start), as.integer(n.cond),
                    as.double(DoF), as.double(cov.chol.sub), start = integer(n.sim.start * n.cond))$start

      dummy <- matrix(dummy, n.sim.start, n.cond, byrow = TRUE)
      dummy.fact <- factor(apply(dummy, 1, paste, collapse = ""))
      start <- dummy[which.max(table(dummy.fact)),]

      cat("Starting partitiion is: ", start, "\n")
      if (model == "Brown-Resnick")
        parts <- .C("gibbsForPartBR", as.integer(k), as.integer(thin), as.integer(burnin),
                as.integer(n.cond), as.integer(start), as.double(cov.sub),
                as.double(sigma2.sub), as.double(cov.chol.sub), as.double(ham.sub),
                as.double(mean1.sub), as.double(y), chain = integer(k * n.cond),
                time1 = double(1))

      if (model == "Schlather")
        parts <- .C("gibbsForPartSC", as.integer(k), as.integer(thin), as.integer(burnin),
                    as.integer(n.cond), as.integer(start), as.double(cov.sub),
                    as.double(y), chain = integer(k * n.cond),
                    time1 = double(1))

      if (model == "extremal-t")
        parts <- .C("gibbsForPartExtt", as.integer(k), as.integer(thin), as.integer(burnin),
                    as.integer(n.cond), as.integer(start), as.double(DoF), as.double(cov.sub),
                    as.double(y), chain = integer(k * n.cond),
                    time1 = double(1))

      timings[1] <- parts$time1
      parts <- parts$chain
      parts <- matrix(parts, n.cond, k)
      cat("    Gibbs partition is: ", parts, "\n")
    }
  }

  else {
    k <- ncol(parts)
    weights <- NULL
  }

  if (do.sim){
    ## Simulation
    if (model == "Brown-Resnick")
      ans <- .C("condsimbrown", as.integer(k), as.integer(n), as.integer(n.cond),
                as.integer(parts), as.double(cov.chol), as.double(sigma2), as.double(ham),
                as.double(mean1), as.double(y.tilde), as.double(all.coord), as.double(range),
                as.double(smooth), as.integer(dim), sim = double(k * n),
                sub.ext.fct = double(k * n), ext.fct = double(k * n), timings = double(2))

    if (model == "Schlather")
      ans <- .C("condsimschlather", as.integer(k), as.integer(n), as.integer(n.cond),
                as.integer(parts), as.double(cov), as.double(y),
                sim = double(k * n), sub.ext.fct = double(k * n),
                ext.fct = double(k * n), timings = double(2))

    if (model == "extremal-t")
      ans <- .C("condsimextt", as.integer(k), as.integer(n), as.integer(n.cond),
                as.integer(parts), as.double(DoF), as.double(cov), as.double(y),
                sim = double(k * n), sub.ext.fct = double(k * n),
                ext.fct = double(k * n), timings = double(2))

    timings[2:3] <- ans$timings
    sub.ext.fct <- matrix(ans$sub.ext.fct, k, n, byrow = TRUE)
    ext.fct <- matrix(ans$ext.fct, k, n, byrow = TRUE)
    ans <- matrix(ans$sim, k, n, byrow = TRUE)

    ## Beware the first n.cond elements correspond to the
    ##conditionning locations!!! If you don't want them then you need
    ##to uncomment the next line...

    ans <- ans[,-(1:n.cond)]
  }

  else
    ans <- sub.ext.fct <- ext.fct <- NULL

  return(list(sim = ans, sub.ext.fct = sub.ext.fct, ext.fct = ext.fct, timings = timings,
              parts = parts))
}
