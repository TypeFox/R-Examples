PRInitialValue2ndLevelMeanKnown <- function(given) {
  # This function makes the initial values needed to run PRIMM.
  # "Kn" means the descriptive second level mean (mean of Beta distribution) is known.
 
  if (all(given$sample.mean == mean(given$sample.mean))) {
    r.ini <- mean(given$prior.mean) / (var(given$sample.mean) + 1)
  } else {
    r.ini <- mean(given$prior.mean) / var(given$sample.mean)
  }

  list(a.ini = -log(r.ini))
}

PRInitialValue2ndLevelMeanUnknown <- function(given) {
  # This function makes the initial values needed to run PRIMM.
  # "Un" means the descriptive second level mean (mean of Beta distribution) is unknown.

  y <- given$sample.mean
  z <- given$z
  n <- given$n
  x.ini <- given$x.ini

  if (identical(x.ini, NA) & !given$intercept) {
    print("In case we do not designate prior mean, at least one beta should be estimated 
           (either intercept=T or at least one covariate is needed)")
    stop()
  } else if (identical(x.ini, NA) & given$intercept) {
    x <- matrix(1, length(y), 1)
    b.ini <- log(mean(y))
  } else if (!identical(x.ini, NA) & given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("z ~ ", paste(xname, collapse = "+")))
    b.ini <- glm(formula, offset = log(n), family = poisson)$coefficients
    x <- as.matrix(cbind(rep(1, length(y)), x.ini))
  } else if (!identical(x.ini,NA) & !given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("z ~ ", paste(xname, collapse = "+"), "- 1"))
    b.ini <- glm(formula, offset = log(n), family = poisson)$coefficients
    x <- x.ini
  }	

  mu0.ini <- mean(exp(x %*% b.ini))

  if (all(y == mean(y))) {
    r.ini <- mu0.ini / (var(y) + 1)
  } else {
    r.ini <- mu0.ini / var(y)
  }

  list(x = x, b.ini = b.ini, a.ini = -log(r.ini))
}
PRAlphaEst2ndLevelMeanKnown <- function(given, ini) {

  z <- given$z
  n <- given$n
  mu0 <- given$prior.mean
  k <- length(n)
  a.ini <- ini$a.ini

  PRDerivAlpha <- function(a) {
    zam <- z + exp(-a) * mu0
    am <- exp(-a) * mu0
    const1 <- ((digamma(zam) - digamma(am) + n / (exp(-a) + n) - log(1 + n * exp(a))) * mu0
               - z / (exp(-a) + n))
    const3 <- (const1 + z * exp(-a) / (exp(-a) + n)^2 + mu0 * n / (n + exp(-a))
              - n * am / (exp(-a) + n)^2 + am * mu0 * (trigamma(zam) - trigamma(am)))
    out <- c(1 - exp(-a) * sum(const1), exp(-a) * sum(const3))
    out
  }

  dif <- 1
  eps <- 0.0001
  while (abs(dif) > eps) {
    out1 <- PRDerivAlpha(a.ini)
    score <- out1[1]
    hessian <- out1[2]
    dif <- score / hessian
    a.ini <- a.ini - dif
  }

  list(a.new = a.ini, a.var = - 1 / hessian)
}

PRAlphaBetaEst2ndLevelMeanUnknown <- function(given, ini) {
  # Alpha and Beta estimation of PRIMM when the descriptive second level mean is unknown

  z <- given$z
  n <- given$n
  x <- ini$x
  k <- length(n)
  a.ini <- ini$a.ini
  b.ini <- ini$b.ini
  m <- ncol(x)

  PRLogLikUn <- function (a, b) { 
    mu0.hat <- exp(x %*% as.matrix(b))
    zam <- z + exp(-a) * mu0.hat
    am <- exp(-a) * mu0.hat
    if (any(c(am, zam) <= 0)) {
      print("The components of lgamma should be positive")
      stop()
    } else {
      sum(dnbinom(z, size = am, prob = exp(-a) / (exp(-a) + n), log=T))
    }
  }

  PRDerivBeta <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    mu0 <- exp(x %*% as.matrix(b))
    zam <- z + exp(-a) * mu0
    am <- exp(-a) * mu0
    if (any(c(am, zam) == 0)) {
      vec[(am == 0 | zam == 0)] <- 0
      tmp <- !(am == 0 | zam == 0)
      vec <- (digamma(zam[tmp]) - digamma(am[tmp]) - log(1 + n[tmp] * exp(a))) * am[tmp]
    } else {
      vec <- (digamma(zam) - digamma(am) - log(1 + n * exp(a))) * am
    }
    if (any(c(am, zam) == 0)) {
      diag[(am == 0 | zam == 0)] <- 0
      tmp <- !(am == 0 | zam == 0)
      diag[tmp] <- (trigamma(zam[tmp]) - trigamma(am[tmp])) * am[tmp]^2 + vec[tmp]
    } else {
      diag <- (trigamma(zam) - trigamma(am)) * am^2 + vec
    }

    out <- cbind(t(x) %*% as.vector(vec), t(x) %*% diag(as.numeric(diag)) %*% x)
    out
  }

  PRDerivBeta2order <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    mu0 <- exp(x %*% as.matrix(b))
    zam <- z + exp(-a) * mu0
    am <- exp(-a) * mu0
    if (any(c(am, zam) == 0)) {
      diag[(am == 0 | zam == 0)] <- 0
      tmp <- !(am == 0 | zam == 0)
      diag[tmp] <- ((trigamma(zam[tmp]) - trigamma(am[tmp])) * am[tmp]^2 
                    + (digamma(zam[tmp]) - digamma(am[tmp]) - log(1 + n[tmp] * exp(a))) * am[tmp])
    } else {
      diag <- ((trigamma(zam) - trigamma(am)) * am^2 
               + (digamma(zam) - digamma(am) - log(1 + n * exp(a))) * am)
    }

    out <- t(x) %*% diag(as.numeric(diag)) %*% x
    out
  }

  BetaHatSubAlpha <- function(a) {
    dif <- 1
    eps <- 0.0001
    n.iter <- 1
    while (any(abs(dif) > eps)) { 
      out <- PRDerivBeta(a, b.ini)
      score <- out[, 1]
      hessian <- out[, 2 : (m + 1)]
      dif <- solve(hessian) %*% score
      b.ini <- b.ini - dif
      n.iter <- n.iter + 1
      if (n.iter > 50) {
        stop()  
        print("Alpha estimate does not converge in Newton-Raphson")
      }
    }
    list(beta.new = b.ini, beta.hessian = PRDerivBeta2order(a, b.ini))
  }

  MarginalPostAlpha <- function(a) {
    b.sub.a <-  BetaHatSubAlpha(a)$beta.new
    a + PRLogLikUn(a, b.sub.a) - 0.5 * log(det(-PRDerivBeta2order(a, b.sub.a)))
  }

  a.temp <- optim(a.ini, MarginalPostAlpha, control = list(fnscale = -1), method= "L-BFGS-B",
                  hessian = TRUE,  lower = -Inf, upper = Inf)
  a.new <- a.temp$par
  a.hess <- a.temp$hessian
  b.temp.result <- BetaHatSubAlpha(a.new)
  b.new <- b.temp.result$beta.new
  b.hessian <- b.temp.result$beta.hessian

  list(a.new = a.new, beta.new = b.new, 
       a.var = -1 / a.hess, beta.var = -solve(b.hessian))
}

PRShrinkageEst <- function(a.res, given) {	
  # This function calculates the shrinkage-related estimates

  a.new <- a.res$a.new
  a.var <- a.res$a.var

  B.hat <- exp(-a.new) / (given$n + exp(-a.new))  # shriankge = B
  inv.info <- 1 / a.var
  var.B.hat <- (B.hat * (1 - B.hat))^2 / ((B.hat * (1 - B.hat)) + inv.info)
  a1.beta <- inv.info / (1 - B.hat)
  a0.beta <- inv.info / B.hat
  B.hat.low <- qbeta((1 - given$confidence.lvl) / 2, a1.beta, a0.beta)
  B.hat.upp <- qbeta((1 + given$confidence.lvl) / 2, a1.beta, a0.beta)
  list(B.hat = B.hat, inv.info = inv.info, var.B.hat = var.B.hat)
}

PRPosteriorEst2ndLevelMeanKnown <- function(B.res, given) {
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is known.

  B.hat <- B.res$B.hat
  var.B.hat <- B.res$var.B.hat
  mu0 <- given$prior.mean
  y <- given$sample.mean
  n <- given$n

  lambda.hat <- y - B.hat * (y - mu0)
  var.lambda.hat <- (lambda.hat - B.hat * y + (var.B.hat + B.hat^2) * y 
                     - (var.B.hat + B.hat^2) * mu0) / n + (y - mu0)^2 * var.B.hat
  u.gamma <- lambda.hat^2 / var.lambda.hat
  v.gamma <- lambda.hat / var.lambda.hat
  lambda.hat.low <- qgamma((1 - given$confidence.lvl) / 2, shape = u.gamma, rate = v.gamma)
  lambda.hat.upp <- qgamma((1 + given$confidence.lvl) / 2, shape = u.gamma, rate = v.gamma)
  list(post.mean = lambda.hat, post.sd = sqrt(var.lambda.hat), 
       post.intv.low = lambda.hat.low, post.intv.upp = lambda.hat.upp, prior.mean = mu0)
}

PRPosteriorEst2ndLevelMeanUnknown <- function(B.res, a.res, ini, given){
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is unknown.

  x <- as.matrix(ini$x)
  n <- given$n
  y <- given$sample.mean
  beta.new <- a.res$beta.new
  beta.var <- a.res$beta.var
  B.hat <- B.res$B.hat
  var.B.hat <- B.res$var.B.hat

  xSx <- diag(x %*% beta.var %*% t(x))
  mu0.hat <- exp(x %*% beta.new + xSx / 2)
  var.mu0 <- mu0.hat^2 * (exp(xSx) - 1)
  lambda.hat <- y - B.hat * (y - mu0.hat)
  var.lambda.hat <- ((lambda.hat - B.hat * y + (var.B.hat + B.hat^2) * y 
                      - (var.B.hat + B.hat^2) * mu0.hat) / n
                     + (var.B.hat + B.hat^2) * (y^2 - 2 * y * mu0.hat + var.mu0 + mu0.hat^2)
                     - (B.hat * (y - mu0.hat))^2)
  u.gamma <- lambda.hat^2 / var.lambda.hat
  v.gamma <- lambda.hat / var.lambda.hat
  lambda.hat.low <- qgamma((1 - given$confidence.lvl) / 2, shape = u.gamma, rate = v.gamma)
  lambda.hat.upp <- qgamma((1 + given$confidence.lvl) / 2, shape = u.gamma, rate = v.gamma)
  list(post.mean = lambda.hat, post.sd = sqrt(var.lambda.hat), 
       post.intv.low = lambda.hat.low, post.intv.upp = lambda.hat.upp, prior.mean = mu0.hat)
}

# main function
pr <- function(z, n, X, prior.mean, intercept = TRUE, confidence.lvl = 0.95) {
  # The main function of PRIMM

  if (missing(X)) { 
    X <- NA
  }

  if (missing(prior.mean)) {
    prior.mean <- NA
  }

  given <- list(z = z, n = n, sample.mean = z/n, x.ini = X, 
                prior.mean = prior.mean, intercept = intercept, confidence.lvl = confidence.lvl)

  if (any(is.na(prior.mean))) {
    ini <- PRInitialValue2ndLevelMeanUnknown(given)
  }else{
    ini <- PRInitialValue2ndLevelMeanKnown(given)
  }

  a.res <- if (missing(prior.mean)) {
             PRAlphaBetaEst2ndLevelMeanUnknown(given, ini)
           } else {
             PRAlphaEst2ndLevelMeanKnown(given, ini)
           }

  B.res <- PRShrinkageEst(a.res, given)

  if (any(is.na(prior.mean))) {
    post.res <- PRPosteriorEst2ndLevelMeanUnknown(B.res, a.res, ini, given)
  }else{
    post.res <- PRPosteriorEst2ndLevelMeanKnown(B.res, given)
  }

  output <- list(sample.mean = given$sample.mean, se = given$n, prior.mean = prior.mean, 
                 shrinkage = B.res$B.hat, sd.shrinkage = sqrt(B.res$var.B.hat), 
                 post.mean = post.res$post.mean, post.sd = post.res$post.sd, 
                 prior.mean.hat = post.res$prior.mean, post.intv.low = post.res$post.intv.low, 
                 post.intv.upp = post.res$post.intv.upp, model = "pr", X = X, 
                 beta.new = a.res$beta.new, beta.var = a.res$beta.var, weight = NA,
                 intercept = intercept, a.new = a.res$a.new, a.var = a.res$a.var, confidence.lvl = confidence.lvl)
  output
}
