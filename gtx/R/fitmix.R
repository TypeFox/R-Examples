## tools for Gaussian mixture models

fitmix.simulate <- function(n, p, mu, sigma) {
  ## simulate n observations from length(p) component mixture with proportions p, means mu and variances sigma^2
  n <- as.integer(n)
  stopifnot(n > 0)
  stopifnot(length(mu) == length(p))
  if (length(sigma) == 1) sigma <- rep(sigma, length(p))
  stopifnot(length(sigma) == length(p))
  m <- sample(1:length(p), n, replace = TRUE, prob = p) # component membership indices
  return(x = mu[m] + rnorm(n)*sigma[m])
}

fitmix.r2 <- function(p, mu, sigma) {
  ## true r2 explained by components in normal mixture, with equal variances for all components
  stopifnot(length(p) == length(mu))
  stopifnot(length(sigma) == 1)
  ve <- sum(p*mu^2)/sum(p) - (sum(p*mu)/sum(p))^2 # Variance explained
  return(r2 = ve/(ve + sigma^2))
}

fitmix.plot <- function(x, p, mu, sigma) {
  ## make nice plot using kernel estimate for observations, plus contribution of mixture
  ## components to mixture density
  stopifnot(length(mu) == length(p))
  if (length(sigma) == 1) sigma <- rep(sigma, length(p))
  stopifnot(length(sigma) == length(p))
  k <- length(p)
  xd <- density(x)
  np <- length(xd$x)
  plot(xd$x, xd$y, type = "n",
       ylim = c(0, max(xd$y*1.5)),
       xlab = "x", ylab = "Density", yaxt = "n")
  rug(x)
  polygon(c(xd$x[1], xd$x, xd$x[np], xd$x[1]),
          c(0, xd$y, 0, 0), col = "grey", lty = 0)
  for (i in 1:k) lines(xd$x, p[i]*dnorm(xd$x, mu[i], sigma[i]), col = rainbow(k)[i], lwd = 3)
  lines(xd$x, rowSums(sapply(1:k, function(i) p[i]*dnorm(xd$x, mu[i], sigma[i]))), lwd = 3)
  legend("topleft", col = c("grey", rainbow(k), "black"), lty = c(0, rep(1, k+1)), lwd = 3,
         pch = c(15, rep(NA, k+1)), pt.cex = 2,
         legend = c("observed", paste("mixture component", 1:k), "mixture total"), bty = "n")
  return(invisible(NULL))
}
  
fitmix1 <- function(x, k, tol = 1e-6, maxit = 100,
                    p = NULL, mu = NULL, sigma = NULL,
                    p.binomial = FALSE, mu.additive = FALSE, sigma.common = FALSE) {
  ## fit mixture of k normal distributions to x using EM
  ## some special options (restricted M steps) for quantitative genetics applications
  x <- as.vector(x)
  x <- x[!is.na(x)] # hard drop NAs 
  stopifnot(length(x) > 1)
  k <- as.integer(k); stopifnot(k > 0)
  tol <- as.double(tol); stopifnot(tol > 0)
  maxit <- as.integer(maxit); stopifnot(maxit > 0)
  if (is.null(p)) p <- rep(1/k, k) # reasonable initial values if not specified
  p <- as.double(p)
  stopifnot(length(p) == k); stopifnot(all(p >= 0)); stopifnot(sum(p) > 0)
  p <- p/sum(p) # normalise if necessary
  if (is.null(mu)) mu <- quantile(x, seq(from = 1/(k+1), to = 1-1/(k+1), length.out = k)) # reasonable initial values if not specified
  mu <- as.double(mu)
  stopifnot(length(mu) == k)
  if (is.null(sigma)) sigma <- rep(sd(x)/3, k) # reasonable initial values if not specified
  sigma <- as.double(sigma)
  stopifnot(length(sigma) == k); stopifnot(all(sigma > 0))
  p.binomial <- as.logical(p.binomial)
  mu.additive <- as.logical(mu.additive)
  sigma.common <- as.logical(sigma.common)
  if (mu.additive & !sigma.common) stop("mu.additive == TRUE requires sigma.common == TRUE")

  ## store progress of loglik by iteration, relative to llnull == loglik for k=1 model
  loglik <- rep(NA, maxit)
  llnull <- -length(x)/2 * (log(2*pi) + log(mean(x^2)-mean(x)^2) + 1)
  
  for (ii in 1:maxit) {
    lij <- sapply(1:k, function(i) return(p[i] * dnorm(x, mu[i], sigma[i])))
    ## likelihood for i-th observation in j-th component 
    loglik[ii] <- sum(log(apply(lij, 1, sum))) - llnull
    if (ii > 1 && loglik[ii] - loglik[ii-1] < tol) break
    pji <- apply(lij, 1, function(pp) return(pp/sum(pp))) # note transposed wrt lij
    ## probability of i-th observation being in j-th component, conditional on p, mu, sigma

    p <- rowSums(pji)/sum(pji) # should have sum(pji)==length(x)
    if (p.binomial) {
      pp <- sum(p*(0:(k-1)))/(k-1)
      p <- choose(k-1,0:(k-1))*pp^(0:(k-1))*(1 - pp)^(k-1:k) # k-1:k == k-(1:k)
    }
    if (!sigma.common) {
      # to ensure correctness, ignore mu.additive argument
      mu <- sapply(1:k, function(j) return(sum(pji[j, ]*x)/sum(pji[j, ]))) # == weighted.mean(x, pji[j, ])
      sigma <- sqrt(sapply(1:k, function(j) sum(pji[j, ]*(x - mu[j])^2)/sum(pji[j, ]))) # estimate for each component
    } else {
      if (!mu.additive) {
        mu <- sapply(1:k, function(j) return(sum(pji[j, ]*x)/sum(pji[j, ]))) # == weighted.mean(x, pji[j, ])
      } else {
        ptmp1 <- rowSums(pji)
        ptmp2 <- sapply(1:k, function(j) return(sum(pji[j, ]*x)))
        stmp <- solve(matrix(c(sum(ptmp1), sum(ptmp1*(0:(k-1))), sum(ptmp1*(0:(k-1))), sum(ptmp1*(0:(k-1))^2)), ncol = 2),
                      matrix(c(sum(ptmp2), sum(ptmp2*(0:(k-1)))), ncol = 1))
        mu <- stmp[1] + (0:(k-1))*stmp[2]
      }
      sigma <- rep(sqrt(sum(sapply(1:k, function(j) sum(pji[j, ]*(x - mu[j])^2)))/sum(pji)), k) # common estimate
    }
  }
  return(list(p = p, mu = mu, sigma = sigma, loglik = loglik[!is.na(loglik)]))
}


fitmix <- function(x, k, tol = 1e-6, maxit = 100, restarts = 20,
                   p.binomial = FALSE, mu.additive = FALSE, sigma.common = FALSE) {
  ## wrapper for fitmix1, calling once with default initial values
  ## then (restarts) times with random initial values
  x <- as.vector(x)
  x <- x[!is.na(x)] # hard drop NAs 
  stopifnot(length(x) > 1)
  k <- as.integer(k); stopifnot(k > 0)
  tol <- as.double(tol); stopifnot(tol > 0)
  maxit <- as.integer(maxit); stopifnot(maxit > 0)
  restarts <- as.integer(restarts); stopifnot(restarts >= 0)
  p.binomial <- as.logical(p.binomial)
  mu.additive <- as.logical(mu.additive)
  sigma.common <- as.logical(sigma.common)
  if (mu.additive & !sigma.common) stop("mu.additive == TRUE requires sigma.common == TRUE")

  afit <- fitmix1(x, k, tol = tol, maxit = maxit,
                  p.binomial = p.binomial, mu.additive = mu.additive, sigma.common = sigma.common)
  p.best <- afit$p
  mu.best <- afit$mu
  sigma.best <- afit$sigma
  loglik <- c(max(afit$loglik), rep(NA, restarts))
  
  if (restarts >= 1) {
    for (run in 1:restarts) {
      ptmp <- rexp(k); ptmp <- ptmp/sum(ptmp) ## Dirichlet(1,1,...,1)
      afit <- fitmix1(x, k, tol = tol, maxit = maxit,
                      p = ptmp, mu = sample(x, 3), sigma = rep(sd(x)/3, 3),
                      p.binomial = p.binomial, mu.additive = mu.additive, sigma.common = sigma.common)
      if (max(afit$loglik) > max(loglik, na.rm = TRUE)) {
        p.best <- afit$p
        mu.best <- afit$mu
        sigma.best <- afit$sigma
      }
      loglik[run + 1] <- max(afit$loglik)
    }
  }
  return(list(p = p.best, mu = mu.best, sigma = sigma.best, loglik = loglik))
}
