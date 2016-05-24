#' @title Simulate Data from DZIM
#' @description Simulate data from a dynamic zero-inflated model. 
#' @param X design matrix.
#' @param w \code{log(w)} is used as an offset variable in the linear predictor.
#' @param omega zero-inflation parameter.
#' @param k dispersion parameter.
#' @param beta regression coefficients.
#' @param phi autoregressive coefficients.
#' @param sigma standard deviation.
#' @param mu0 mean vector of initial state.
#' @param Sigma0 covariance matrix of initial state.
#' @seealso 
#' \code{\link{dzim}},
#' \code{\link{dzim.fit}},
#' \code{\link{dzim.filter}}, 
#' \code{\link{dzim.smooth}},
#' \code{\link{dzim.control}},
#' \code{\link{dzim.plot}}
#' @keywords regression
#' @export dzim.sim
dzim.sim <-  function(X, w, omega, k, beta, phi, sigma, mu0, Sigma0) {
  n <- NROW(X)
  p <- length(phi)
  s0 <- mvrnorm(1, mu0, Sigma0)
  s <- matrix(NA, n, p)
  u <- rep(NA, n) 
  v <- rep(NA, n) 
  y <- rep(NA, n)
  Phi <- suppressWarnings(rbind(phi, cbind(diag(1, p - 1), 0))) 
  s[1, ] <- Phi %*% s0 +  c(rnorm(1, 0, sigma), rep(0, p - 1))
  u[1] <- rbinom(1, 1, omega)
  v[1] <- ifelse(k == Inf, 1, rgamma(1, shape = k, scale = 1 / k))
  lambda <- w[1] * exp(sum(X[1, ] * beta) + s[1, 1])                        
  y[1] <- rpois(1, (1 - u[1]) * v[1] * lambda)
  for(t in 2:n) {
    s[t, ] <- Phi %*% s[t - 1, ] + c(rnorm(1, 0, sigma), rep(0, p - 1))
    u[t] <- rbinom(1, 1, omega)
    v[t] <- ifelse(k == Inf, 1, rgamma(1, shape = k, scale = 1 / k))
    lambda <- w[t] * exp(sum(X[t, ] * beta) + s[t, 1])
    y[t] <- rpois(1, (1 - u[t]) * v[t] * lambda)
  }
  list(s0 = s0, s = s, u = u, v = v, y = y)
}

#' @title Particle Filtering for DZIM
#' @description Function to implement the particle filtering method proposed by Gordsill et al. (1993). 
#' @param y response variable.
#' @param X design matrix.
#' @param w \code{log(w)} is used as an offset variable in the linear predictor.
#' @param para model parameters.
#' @param control control arguments.
#' @seealso
#' \code{\link{dzim}},
#' \code{\link{dzim.fit}}, 
#' \code{\link{dzim.smooth}},
#' \code{\link{dzim.control}},
#' \code{\link{dzim.sim}},
#' \code{\link{dzim.plot}}
#' @references
#' Gordon, N. J., Salmond, D. J., and Smith, A. F. M. (1993). Novel approach to nonlinear/non-Gaussian 
#' Bayesian state estimation. \emph{IEEE Proceedings}, \bold{140}, 107-113. 
#' @keywords regression
#' @export dzim.filter
dzim.filter <- function(y, X, w, para, control) {
  n <- NROW(X)
  pX <- NCOL(X)
  p <- control$order
  N <- control$N
  R <- control$R
  mu0 <- control$mu0
  Sigma0 <- control$Sigma0
  para <- switch(control$dist, "poisson" = c(0, Inf, para), 
    "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
  omega <- para[1]
  k <- para[2]
  beta <- para[2 + 1:pX]
  phi <- para[2 + pX + 1:p]
  sigma <- para[length(para)]
  sp <- array(NA, dim = c(N, n, p)) 
  up <- matrix(NA, N, n)
  vp <- matrix(NA, N, n)
  qf <- matrix(NA, N, n) 
  sf <- array(NA, dim = c(N, n, p)) 
  uf <- matrix(NA, N, n) 
  vf <- matrix(NA, N, n)
  Phi <- suppressWarnings(rbind(phi, cbind(diag(1, p - 1), 0)))  
  qf0 <- rep(1, N)
  sf0 <- mvrnorm(N, mu0, Sigma0)
  for(i in 1:N) {
    sp[i, 1, ] <- Phi %*% sf0[i, ] + c(rnorm(1, 0, sigma), rep(0, p - 1))
  }
  up[, 1] <- rbinom(N, 1, omega)
  vp[, 1] <- ifelse(rep(k == Inf, N), rep(1, N), rgamma(N, shape = k, scale = 1 / k))
  lambda <- w[1] * exp(sum(X[1, ] * beta) + sp[, 1, 1])
  qf[, 1] <- dpois(y[1], (1 - up[, 1]) * vp[, 1] * lambda)
  index <- suppressWarnings(sample(1:N,
    size = N, replace = TRUE, prob = qf[, 1]))
  sf[, 1, ] <- sp[index, 1, ]
  uf[, 1] <- up[index, 1]
  vf[, 1] <- vp[index, 1]
  for(t in 2:n) {
    for(i in 1:N) {
      sp[i, t, ] <- Phi %*% sf[i, t - 1, ] + c(rnorm(1, 0, sigma), rep(0, p - 1))    
    }
    up[, t] <- rbinom(N, 1, omega)
    vp[, t] <- ifelse(rep(k == Inf, N), rep(1, N), rgamma(N, shape = k, scale = 1 / k))
    lambda <- w[t] * exp(sum(X[t, ] * beta) + sp[, t, 1])
    qf[, t] <- dpois(y[t], (1 - up[, t]) * vp[, t] * lambda)
    index <- suppressWarnings(sample(1:N,
      size = N, replace = TRUE, prob = qf[, t]))
    sf[, t, ] <- sp[index, t, ]
    uf[, t] <- up[index, t]
    vf[, t] <- vp[index, t]
  }
  loglik <- sum(log(colMeans(qf)))
  list(qf0 = qf0, sf0 = sf0, sp = sp, up = up, vp = vp, 
    qf = qf, sf = sf, uf = uf, vf = vf, loglik = loglik)
}

#' @title Particle Smoothing for DZIM
#' @description Function to implement the particle smoothing method proposed by Gordsill et al. (2004).
#' @param y response variable.
#' @param X design matrix.
#' @param w \code{log(w)} is used as an offset variable in the linear predictor.
#' @param para model parameters.
#' @param control control arguments.
#' @seealso
#' \code{\link{dzim}},
#' \code{\link{dzim.fit}}, 
#' \code{\link{dzim.filter}},
#' \code{\link{dzim.control}}, 
#' \code{\link{dzim.sim}}, 
#' \code{\link{dzim.plot}}
#' @references 
#' Gordsill, S. J., Doucet, A., and West, M. (2004). Monte Carlo smoothing for nonlinear time series. 
#' \emph{Journal of the American Statistical Association}, \bold{99}, 156-168. 
#' @keywords regression
#' @export dzim.smooth
dzim.smooth <- function(y, X, w, para, control) {
  n <- NROW(X)
  pX <- NCOL(X)
  p <- control$order
  N <- control$N
  R <- control$R
  mu0 <- control$mu0
  Sigma0 <- control$Sigma0
  pf <- dzim.filter(y, X, w, para, control)
  para <- switch(control$dist, "poisson" = c(0, Inf, para), 
    "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
  omega <- para[1]
  k <- para[2]
  beta <- para[2 + 1:pX]
  phi <- para[2 + pX + 1:p]
  sigma <- para[length(para)] 
  qs <- matrix(NA, N, n) 
  ss <- array(NA, dim = c(R, n, p)) 
  us <- matrix(NA, R, n) 
  vs <- matrix(NA, R, n)
  qs0 <- rep(NA, N)
  ss0 <- matrix(NA, R, p) 
  for(r in 1:R) {
    qs[, n] <- pf$qf[, n]
    index <- sample(1:N, size = 1, prob = qs[, n])
    ss[r, n, ] <- pf$sf[index, n, ]
    us[r, n] <- pf$uf[index, n]
    vs[r, n] <- pf$vf[index, n]
    for(t in (n - 1):1) {
      qs[, t] <- pf$qf[, t] *
        dnorm(ss[r, t + 1, 1], as.matrix(pf$sf[, t, ]) %*% phi, sigma) *
        dbinom(us[r, t + 1], 1, omega) *
        ifelse(k == Inf, 1, dgamma(vs[r, t + 1], shape = k, scale = 1 / k))
      index <- sample(1:N, size = 1, prob = qs[, t])
      ss[r, t, ] <- pf$sf[index, t, ]
      us[r, t] <- pf$uf[index, t]
      vs[r, t] <- pf$vf[index, t]    
    }
    qs0 <- pf$qf0 *
        dnorm(ss[r, 1, 1], as.matrix(pf$sf0) %*% phi, sigma) *
        dbinom(us[r, 1], 1, omega) *
        ifelse(k == Inf, 1, dgamma(vs[r, 1], shape = k, scale = 1 / k))
    index <- sample(1:N, size = 1, prob = qs0)    
    ss0[r, ] <- pf$sf0[index, ]  
  }
  list(pf = pf, qs = qs, ss = ss, us = us, vs = vs, qs0 = qs0, ss0 = ss0)  
}

#' @title Auxiliary for Controlling DZIM Fitting
#' @description Auxiliary function for \code{\link{dzim}} fitting. Typically only used internally by 
#' \code{\link{dzim.fit}}, but may be used to construct a control argument for either function.
#' @param dist count model family
#' @param trace logical; if TRUE, display iteration history.
#' @param start initial parameter values.
#' @param order autoregressive order.
#' @param mu0 mean vector for initial state.
#' @param Sigma0 covariance matrix for initial state.
#' @param N number of particiles in particle filtering.
#' @param R number of replications in particle smoothing.
#' @param niter number of iterations.
#' @note The default values of \code{N}, \code{R}, and \code{niter} are chosen based on our experience.
#' In some cases, \code{N} = 500, \code{R} = 500, and \code{niter} = 200 might be sufficient. 
#' The \code{\link{dzim.plot}} function should always be used for convergence diagnostics. 
#' @seealso
#' \code{\link{dzim}},
#' \code{\link{dzim.fit}}, 
#' \code{\link{dzim.filter}}, 
#' \code{\link{dzim.smooth}},
#' \code{\link{dzim.sim}},
#' \code{\link{dzim.plot}}
#' @keywords regression
#' @export dzim.control
dzim.control <- function(dist = c("poisson", "nb", "zip", "zinb"),
  trace = FALSE, start = NULL, order = 1, 
  mu0 = rep(0, order), Sigma0 = diag(1, order), 
  N = 1000, R = 1000, niter = 500) {
  dist <- match.arg(dist)
  list(dist = dist, trace = trace, start = start,
    order = order, mu0 = mu0, Sigma0 = Sigma0,
    N = N, R = R, niter = niter)
}

#' @title Fitter Function for Dynamic Zero-Inflated Models
#' @description \code{\link{dzim.fit}} is the basic computing engine called by \code{\link{dzim}} used to fit 
#' dynamic zero-inflated models. This should usually \emph{not} be used directly unless by experienced users. 
#' @param y response variable.
#' @param X design matrix.
#' @param offset offset variable.
#' @param control control arguments.
#' @param ... additional arguments.
#' @seealso
#' \code{\link{dzim}},
#' \code{\link{dzim.control}}, 
#' \code{\link{dzim.filter}}, 
#' \code{\link{dzim.smooth}},
#' \code{\link{dzim.sim}},
#' \code{\link{dzim.plot}}
#' @keywords regression
#' @export dzim.fit
dzim.fit <- function(y, X, offset = rep(0, n), control = dzim.control(...), ...) {
  deriv <- function(para) {
    R <- control$R
    para <- switch(control$dist, "poisson" = c(0, Inf, para), 
      "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
    omega <- para[1]
    k <- para[2]
    beta <- para[2 + 1:pX]
    phi <- para[2 + pX + 1:p]
    sigma <- para[length(para)]
    m <- length(para)
    score <- rep(0, m)
    info.com <- matrix(0, m, m)
    info.mis <- matrix(0, m, m)
    info.obs <- matrix(0, m, m)
    U <- array(NA, dim = c(m, m, R))
    V <- array(NA, dim = c(m, m, R))
    W <- array(NA, dim = c(m, m, R))
    for(r in 1:R) {
      temp <- rep(0, m)
      v <- rep(0, m)
      K <- matrix(0, m, m)
      M <- matrix(0, m, m)  
      v[1] <- ifelse(omega == 0, 0, 
        ps$us[r, 1] / omega - (1 - ps$us[r, 1]) / (1 - omega))
      v[2] <- ifelse(k == Inf, 0, 
        1 + log(k) - digamma(k) + log(ps$vs[r, 1]) - ps$vs[r, 1])
      v[2 + 1:pX] <- (1 - ps$us[r, 1]) * (y[1] - 
        ps$vs[r, 1] * w[1] * exp(sum(X[1, ] * beta) + ps$ss[r, 1, 1])) * X[1, ] 
      v[2 + pX + 1:p] <- (1 / sigma^2) * 
        (ps$ss[r, 1, 1] - sum(phi * ps$ss0[r, ])) * ps$ss0[r, ]  
      v[m] <- (ps$ss[r, 1, 1] - 
        sum(phi * ps$ss0[r, ]))^2 / sigma^3 - 1 / sigma            
      K <- v %*% t(v)
      M[1, 1] <- ifelse(omega == 0, 0, 
        ps$us[r, 1] / omega^2 + (1 - ps$us[r, 1]) / (1 - omega)^2)
      M[2, 2] <- ifelse(k == Inf, 0, 
        trigamma(k) - 1 / k)
      M[2 + 1:pX, 2 + 1:pX] <- (1 - ps$us[r, 1]) * ps$vs[r, 1] * w[1] * 
        exp(sum(X[1, ] * beta) + ps$ss[r, 1, 1]) * X[1, ] %*% t(X[1, ])
      M[2 + pX + 1:p, 2 + pX + 1:p] <- (1 / sigma^2) * 
        ps$ss0[r, ] %*% t(ps$ss0[r, ])
      M[m, m] <- (3 / sigma^4) *  
        (ps$ss[r, 1, 1] - sum(phi * ps$ss0[r, ]))^2 - 1 / sigma^2
      M[2 + pX + 1:p, m] <- (2 / sigma^3) *
        (ps$ss[r, 1, 1] - sum(phi * ps$ss0[r, ])) * ps$ss0[r, ]               
      for(t in 2:n) {
        temp[1] <- ifelse(omega == 0, 0, 
          ps$us[r, t] / omega - (1 - ps$us[r, t]) / (1 - omega))
        temp[2] <- ifelse(k == Inf, 0, 
          1 + log(k) - digamma(k) + log(ps$vs[r, t]) - ps$vs[r, t])
        temp[2 + 1:pX] <- (1 - ps$us[r, t]) * (y[t] - 
          ps$vs[r, t] * w[t] * exp(sum(X[t, ] * beta) + ps$ss[r, t, 1])) * X[t, ] 
        temp[2 + pX + 1:p] <- (1 / sigma^2) * 
          (ps$ss[r, t, 1] - sum(phi * ps$ss[r, t - 1, ])) * ps$ss[r, t - 1, ] 
        temp[m] <- (ps$ss[r, t, 1] - 
          sum(phi * ps$ss[r, t - 1, ]))^2 / sigma^3 - 1 / sigma 
        v[1] <- v[1] + temp[1]
        v[2] <- v[2] + temp[2]
        v[2 + 1:pX] <- v[2 + 1:pX] + temp[2 + 1:pX]
        v[2 + pX + 1:p] <- v[2 + pX + 1:p] + temp[2 + pX + 1:p]
        v[m] <- v[m] + temp[m]
        K <- K + temp %*% t(temp)   
        M[1, 1] <- M[1, 1] + ifelse(omega == 0, 0, 
          ps$us[r, t] / omega^2 + (1 - ps$us[r, t]) / (1 - omega)^2)
        M[2, 2] <- M[2, 2] + ifelse(k == Inf, 0, 
          trigamma(k) - 1 / k)
        M[2 + 1:pX, 2 + 1:pX] <- M[2 + 1:pX, 2 + 1:pX] + 
          (1 - ps$us[r, t]) * ps$vs[r, t] * w[t] * 
          exp(sum(X[t, ] * beta) + ps$ss[r, t, 1]) * X[t, ] %*% t(X[t, ])
        M[2 + pX + 1:p, 2 + pX + 1:p] <- M[2 + pX + 1:p, 2 + pX + 1:p] + 
          (1 / sigma^2) * ps$ss[r, t - 1, ] %*% t(ps$ss[r, t - 1, ])
        M[m, m] <- M[m, m] + (3 / sigma^4) *  
          (ps$ss[r, t, 1] - sum(phi * ps$ss[r, t - 1, ]))^2 - 1 / sigma^2
        M[2 + pX + 1:p, m] <- M[2 + pX + 1:p, m] + (2 / sigma^3) *
          (ps$ss[r, t, 1] - sum(phi * ps$ss[r, t - 1, ])) * ps$ss[r, t - 1, ]       
      }
      M[m, 2 + pX + 1:p] <- t(M[2 + pX + 1:p, m])
      score <- score + v / R
      U[, , r] <- M
      V[, , r] <- v %*% t(v)
      W[, , r] <- K
    }
    gradient <- switch(control$dist, "poisson" = score[-(1:2)], 
      "nb" = score[-1], "zip" = score[-2], "zinb" = score)
    J <- apply(W, c(1, 2), mean)
    J <- switch(control$dist, "poisson" = J[-(1:2), -(1:2)], 
      "nb" = J[-1, -1], "zip" = J[-2, -2], "zinb" = J)
    J <- (J + t(J)) / 2
    info.com <- apply(U, c(1, 2), mean)
    info.mis <- apply(V, c(1, 2), mean) - score %*% t(score)
    for(prop in seq(0, 1, 0.01)) {
      info <- info.com - (1 - prop) * info.mis
      info <- switch(control$dist, "poisson" = info[-(1:2), -(1:2)], 
        "nb" = info[-1, -1], "zip" = info[-2, -2], "zinb" = info)
      info <- (info + t(info)) / 2 
      if(all(eigen(info)$values > 0)) {
        break
      }     
    }
    list(gradient = gradient, J = J, info = info)
  }  
  ar.fit <- function(ss, ss0) {
    ss.dim <- dim(ss)
    R <- ss.dim[1]
    n <- ss.dim[2]
    p <- ss.dim[3]
    A <- array(0, dim = c(R, p, p))
    b <- matrix(0, R, p)
    for(r in 1:R) {
        A[r, , ] <- ss0[r, ] %*% t(ss0[r, ])
        b[r, ] <- ss[r, 1, 1] * ss0[r, ]
      for(t in 2:n) {
        A[r, , ] <- A[r, , ] + ss[r, t - 1, ] %*% t(ss[r, t - 1, ])
        b[r, ] <- b[r, ] + ss[r, t, 1] * ss[r, t - 1, ]
      }
    }
    A <- apply(A, c(2, 3), mean)
    b <- colMeans(b)
    c <- colMeans(ss[, , 1]^2)
    phi <- solve(A) %*% b
    sigma <- sqrt((sum(c) - t(b) %*% solve(A) %*% b) / n)
    list(phi = c(phi), sigma = c(sigma))
  }
  nb.fit <- function(vs) {                     
    h <- function(p) {
      k <- p / (1 - p)
      e <- colMeans(vs)
      f <- colMeans(log(vs))
      abs(sum(1 + log(k) - digamma(k) + f - e))
    }
    p <- optimize(h, c(0, 1))$minimum
    p / (1 - p)
  }
  em <- function(para) {
    para <- switch(control$dist, "poisson" = c(0, Inf, para), 
      "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
    omega <- para[1]
    k <- para[2]
    beta <- para[2 + 1:pX]
    phi <- para[2 + pX + 1:p]
    sigma <- para[length(para)]
    d <- colMeans(ps$us)
    g <- colMeans((1 - ps$us) * ps$vs * exp(ps$ss[, , 1]))  
    omega <- ifelse(sum(control$dist == c("poisson", "nb")), 0, mean(d))
    k <- ifelse(sum(control$dist == c("poisson", "zip")), Inf, nb.fit(ps$vs))
    beta <- glm.fit(X, y, weights = 1 - d, 
     offset = ifelse(1 - d, log(g * w / (1 - d)), 0),
     family = poisson(), start = beta)$coef
    ar <- ar.fit(ps$ss, ps$ss0)  
    phi <- ar$phi
    sigma <- ar$sigma
    switch(control$dist,
      "poisson" = c(beta, phi, sigma),
      "nb" = c(k, beta, phi, sigma),
      "zip" = c(omega, beta, phi, sigma),
      "zinb" = c(omega, k, beta, phi, sigma))
  }
  n <- NROW(X)
  pX <- NCOL(X)
  w <- exp(offset)
  p <- control$order
  iter <- 0
  if(is.null(control$start)) {
    if(control$dist == "poisson") {
      fit1 <- glm(y ~ 0 + X + offset(offset), family = poisson)
      omega.old <- 0
      k.old <- Inf
      beta.old <- fit1$coef
    } else if(control$dist == "nb") {
      fit2 <- glm.nb(y ~ 0 + X + offset(offset))
      omega.old <- 0
      k.old <- fit2$theta
      beta.old <- fit2$coef
    } else if(control$dist == "zip") {
      fit3 <- zim(y ~ 0 + X + offset(offset) | 1, dist = "zip")
      omega.old <- plogis(fit3$para[pX + 1])
      k.old <- Inf
      beta.old <- fit3$para[1:pX]
    } else if(control$dist == "zinb") {
      fit4 <- zim(y ~ 0 + X + offset(offset) | 1, dist = "zinb")
      omega.old <- plogis(fit4$para[1 + pX + 1])  
      k.old <- exp(fit4$para[1])
      beta.old <- fit4$para[1 + 1:pX]    
    }
    phi.old <- rep(0, control$order)
    sigma.old <- 1    
  }
  para.old <- switch(control$dist,
    "poisson" = c(beta.old, phi.old, sigma.old),
    "nb" = c(k.old, beta.old, phi.old, sigma.old),
    "zip" = c(omega.old, beta.old, phi.old, sigma.old),
    "zinb" = c(omega.old, k.old, beta.old, phi.old, sigma.old))
  para.trace <- NULL
  loglik.trace <- NULL
  for(iter in 1:control$niter) {
    ps <- dzim.smooth(y, X, w, para.old, control)
    para.new <- em(para.old)
    para.old <- para.new
    loglik.old <- ps$pf$loglik 
    para.trace <- rbind(para.trace, para.old)
    loglik.trace <- c(loglik.trace, loglik.old) 
    if(control$trace == TRUE) {
      cat("iter =", iter, "\t loglik =", loglik.old, "\n")
    }
  }
  para <- para.new
  para.names <- c("omega", "k", paste("beta", as.character(1:pX - 1), sep = ""),
    paste("ar", as.character(1:control$order), sep = ""), "sigma")
  para.names <- switch(control$dist, "poisson" = para.names[-(1:2)], 
    "nb" = para.names[-1], "zip" = para.names[-2], "zinb" = para.names)
  colnames(para.trace) <- para.names
  names(para) <- para.names
  control$N <- 10 * control$N
  control$R <- 10 * control$R  
  ps <- dzim.smooth(y, X, w, para.new, control)
  deriv <- deriv(para.new)
  para.tmp <- switch(control$dist, "poisson" = c(0, Inf, para), 
    "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
  omega <- para.tmp[1]
  k <- para.tmp[2]
  beta <- para.tmp[2 + 1:pX]
  phi <- para.tmp[2 + pX + 1:p]
  sigma <- para.tmp[length(para)]
  lambda <- exp(offset + X %*% beta) * colMeans(exp(ps$ss[, , 1]))
  mu <- (1 - omega) * lambda
  names(deriv$gradient) <- names(para)
  rownames(deriv$J) <- names(para)
  colnames(deriv$J) <- names(para)
  rownames(deriv$info) <- names(para)
  colnames(deriv$info) <- names(para)
  list(omega = omega, k = k, beta = beta, phi = phi, sigma = sigma, 
    lambda = as.vector(lambda), mu = as.vector(mu),
    para = para, ps = ps, deriv = deriv,
    para.trace = ts(para.trace), loglik.trace = ts(loglik.trace))
}                               

#' @title Fitting Dynamic Zero-Inflated Models
#' @description \code{dzim} is used to fit dynamic zero-inflated models. 
#' @param formula an objective of class "\code{\link{formula}}".
#' @param data an optional dataframe, list or environment containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s.
#' @param weights an optional vector of 'prior weights' to be used in the fitting process.
#' @param offset this can be used to specify a priori known component to be included in the linear predictor during fitting.
#' @param control control arguments from \code{\link{dzim.control}}
#' @param ... additional arguments
#' @seealso
#' \code{\link{dzim.fit}},
#' \code{\link{dzim.filter}}, 
#' \code{\link{dzim.smooth}},
#' \code{\link{dzim.control}}, 
#' \code{\link{dzim.sim}},
#' \code{\link{dzim.plot}}
#' @keywords regression
#' @export dzim
dzim <- function(formula, data, subset, na.action, weights = 1, offset = 0, 
  control = dzim.control(...), ...) {
  call <- match.call()
  if(missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action",
    "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf[[1]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE 
  mf$formula <- formula
  mf <- eval(mf, parent.frame())    
  y <- round(model.response(mf, "numeric"))
  X <- model.matrix(terms(formula, data = data), mf)
  n <- NROW(X)
  weights <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  if(is.null(weights)) {
    weights <- rep(1, n)
  }
  if(is.null(offset)) {
    offset <- rep(0, n)
  }    
  fit <- dzim.fit(y, X, offset = offset, control = control)
  fit$call <- call
  fit$control <- control
  fit$na.action <- attr(mf, "na.action")
  fit$y <- y
  fit$X <- X  
  fit$aic <- (-2) * fit$ps$pf$loglik + 2 * length(fit$para)  
  fit$bic <- (-2) * fit$ps$pf$loglik + log(n) * length(fit$para) 
  if(is.numeric(try(solve(fit$deriv$info)))) {
    fit$se <- sqrt(diag(solve(fit$deriv$info)))
    if(any(is.na(fit$se))) {
      fit$se <- rep(NA, length(fit$para))    
    }
    fit$tic <- (-2) * fit$ps$pf$loglik + 
      2 * sum(diag(fit$deriv$J %*% solve(fit$deriv$info)))    
  } else {
    fit$se <- rep(NA, length(fit$para)) 
    fit$tic <- NA
  }
  names(fit$se) <- names(fit$para)      
  class(fit) <- "dzim"
  fit
}

#' @export
print.dzim <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  pX <- NCOL(x$X)
  z.value <- x$para / x$se 
  z.prob <- pvalue(z.value)
  coef <- data.frame(x$para, x$se, z.value, z.prob)  
  coefX <- switch(x$control$dist, 
    "poisson" = coef[1:pX, ], "nb" = coef[1 + 1:pX, ], 
    "zip" = coef[1 + 1:pX, ], "zinb" = coef[2 + 1:pX, ])
  coefAR<- switch(x$control$dist, 
    "poisson" = coef[-(1:pX), ], "nb" = coef[-(1:(1 + pX)), ], 
    "zip" = coef[-(1:(1 + pX)), ], "zinb" = coef[-(1:(2 + pX)), ])
  coefAR <- coefAR[1:x$control$order, ]
  colnames(coefX) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)") 
  colnames(coefAR) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(coefX) <- colnames(x$X)
  rownames(coefAR) <- paste("ar", 1:x$control$order, sep = "")  
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if(x$control$dist == "nb") {                   
    cat(paste("(Dispersion parameter for negative binomial taken to be ",
      round(x$para[1], 4), ")", sep = ""), "\n\n")
  } else if(x$control$dist == "zip") {
    cat(paste("(Zero-inflation parameter taken to be ",
      round(x$para[1], 4), ")", sep = ""), "\n\n")
  } else if(x$control$dist == "zinb") {
    cat(paste("(Zero-inflation parameter taken to be ",
      round(x$para[1], 4), ")", sep = ""), "\n\n")
    cat(paste("(Dispersion parameter for negative binomial taken to be ",
      round(x$para[2], 4), ")", sep = ""), "\n\n")
  } 
  cat("Coefficients (log-linear): \n")     
  printCoefmat(coefX, signif.legend = FALSE)
  cat("\n")
  cat("Coefficients (autoregressive): \n")
  printCoefmat(coefAR, signif.legend = FALSE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n") 
    cat(paste("(Standard deviation parameter taken to be ",
      round(x$para[length(x$para)], 4), ")", sep = ""), "\n\n")  
  cat("Criteria for assessing goodness of fit \n") 
  cat("loglik:", x$ps$pf$loglik, "\n")
  cat("aic:", x$aic, "\n")
  cat("bic:", x$bic, "\n") 
  cat("tic:", x$tic, "\n") 
  cat("\n")
  invisible(x)
}

#' @title Trace Plots from DZIM
#' @description Function to display trace plots from a dynamic zero-inflated model.
#' @param object objective from \code{\link{dzim}} or \code{\link{dzim.fit}}.
#' @param k.inv logical; indicating whether an inverse transformation is needed for the dispersion parameter.
#' @param sigma.sq logical; indicating whether a square transformation is needed for the standard deviation parameter.
#' @param ... additional arguments.
#' @export dzim.plot
dzim.plot <- function(object, k.inv = FALSE, sigma.sq = FALSE, ...) {
  pX <- NCOL(object$X)
  control <- object$control
  para.trace <- object$para.trace
  I <- NROW(para.trace)
  J <- NCOL(para.trace)
  if(k.inv == TRUE) {
    if(control$dist == "nb")   para.trace[, 1] <- 1 / para.trace[, 1]
    if(control$dist == "zinb") para.trace[, 2] <- 1 / para.trace[, 2]
  }
  if(sigma.sq == TRUE) {
    para.trace[, NCOL(para.trace)] <- para.trace[, NCOL(para.trace)]^2
  }
  scale <- apply(para.trace[round(I/2):I, ], 2, sd)
  scale <- ifelse(scale == 0, 1, scale)
  para.trace.plot <- (para.trace -  matrix(rep(para.trace[1, ], I), nrow = I, byrow = TRUE)) %*% solve(diag(scale))
  col <- c(rep("green", pX), rep("blue", control$order), "purple")
  lty <- c(1:pX, 1:control$order, 1)
  col <- switch(control$dist, "poisson" = col, 
    "nb" = c("red", col), "zip" = c("gold", col), "zinb" = c("gold", "red", col))
  lty <- switch(control$dist, "poisson" = lty, 
    "nb" = c(1, lty), "zip" = c(1, lty), "zinb" = c(1, 1, lty))  
  plot.ts(para.trace.plot, plot.type = "single", lty = lty, col = col, yaxt = "n", xlab = "Iteration", ylab = " ", ...)
  points(1, 0, pch = 20)
}

