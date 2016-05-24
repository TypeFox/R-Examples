
#'@title Fish stable isotope dataset.
#'@description A dataset containing values for three stable isotopes measured
#'in the muscle tissue of four species of arctic fish. For use in examples.
#'@details This dataset contains delta^{15}N, delta^{13}C, and delta^{34}S
#'values for the following fish species:
#'\itemize{
#'\item ARCS - Arctic Cisco (\emph{Coregonus autumnalis}), n = 69
#'\item BDWF - Broad Whitefish (\emph{Coregonus nasus}), n = 71
#'\item LKWF - Lake Whitefish (\emph{Coregonus clupeaformis}), n = 67
#'\item LSCS - Least Cisco (\emph{Coregonus sardinella}), n = 70
#'}
#'Fish were collected between 2007 and 2008 from an estuarine area of the
#'Beaufort Sea, North and West of the Mackenzie Delta at Phillips
#'Bay, Yukon Territory, Canada (69.28 N, 138.49 W).
#'@docType data
#'@format A data frame with 278 rows (observations) and 4 columns (species, delta^{15}N, delta^{13}C, and delta^{34}S).
#'@examples
#'data(fish)
#'aggregate(fish[2:4], fish[1], mean)
#'@name fish
NULL

#################################################

#'@title Point coordinates for a 2-D ellipse.
#'@description Calculates coordinates of points for plotting a 2-dimensional ellipse based on
#'user-defined parameters. Can be used for exploratory data analysis to produce ellipses 
#'at a given niche region size (e.g., \eqn{\alpha = 95\%}).
#'@details This function provides the coordinates needed to plot a 2-dimensional ellipse
#'based on user-defined parameters, such that \code{X = c(x,y)} satisfies the equation
#'\deqn{(X-\mu)' V^{-1} (X-\mu) = C,}
#'where \eqn{C=\code{qchisq(alpha, df = 2)}}.
#'@param mu centre of ellipse. A vector of length 2.
#'@param V scale of ellipse. A 2x2 matrix. See Details.
#'@param alpha niche region size. See Details.
#'@param n number of points to return for plotting.
#'@return Returns a matrix of coordinates \code{cbind(x,y)} to plot a 2-dimensional ellipse.
#'@seealso \code{\link{niche.plot}}
#'@examples
#'mu <- rnorm(2)
#'V <- crossprod(matrix(rnorm(4), 2, 2))
#'ell.pts <- ellipse(mu = mu, V = V, alpha = .9, n = 100)
#'plot(ell.pts, col = rainbow(110)[1:100], type = "o")
#'points(mu[1], mu[2], pch = "+")
#'@export
ellipse <- function(mu, V, alpha = .95, n = 100) {
  tmp <- eigen(V)
  hlen <- sqrt(qchisq(alpha, df = 2)*tmp$val)
  theta <- atan2(tmp$vec[2,1], tmp$vec[1,1])
  t <- seq(0, 2*pi, len = n+1)
  x <- hlen[1] * cos(t)
  y <- hlen[2] * sin(t)
  alpha <- atan2(y, x)
  rad <- sqrt(x^2 + y^2)
  cbind(x = rad * cos(alpha + theta) + mu[1],
        y = rad * sin(alpha + theta) + mu[2])
}

#################################################

#'@title Random draws from a Wishart (or Inverse-Wishart) distribution.
#'@description Generates a random samples from a Wishart distribution defined as
#'\eqn{W(\Psi, \nu)}, or an Inverse-Wishart distribution defined as \eqn{W^{-1}(\Psi, \nu)}. 
#'@details Setting \code{inv = TRUE} replaces \eqn{\Psi} by \eqn{Psi^{-1}} and inverts the output random matrices,
#'such that they are being generated from an Inverse-Wishart \eqn{W^{-1}(\Psi, \nu)} distribution.
#'@param n number of samples to draw.
#'@param Psi scale matrix.
#'@param nu degrees of freedom.
#'@param inv logical. Setting \code{inv = TRUE} returns random matrices from an Inverse-Wishart
#'distribution. See Details.
#'@seealso \code{\link{rniw}}
#'@return Returns an array of Wishart (or Inverse-Wishart) draws of size \code{c(nrow(Psi),ncol(Psi),n)}.
#'@examples
#'d <- 4 # number of dimensions
#'nu <- 7 # degrees of freedom
#'Psi <- crossprod(matrix(rnorm(d^2), d, d)) # scale matrix
#'n <- 1e4
#'
#'Sigma <- rwish(n, Psi, nu)
#'
#'# for any vector a, X = (a' Sigma a) has a const * chi^2 distribution
#'a <- rnorm(d)
#'X <- apply(Sigma, 3, function(S) crossprod(a, S %*% a))
#'const <- a %*% Psi %*% a
#'
#'hist(X, breaks = 100, freq = FALSE,
#'     main = parse(text = "\"Histogram of \"*X==a*minute*Sigma*a"),
#'     xlab = parse(text = "X==a*minute*Sigma*a"))
#'curve(dchisq(x/const, df = nu)/const,
#'      from = min(X), to = max(X), col = "red", add = TRUE)
#'@export
rwish <- function(n, Psi, nu, inv = FALSE) {
  if(inv) Psi <- solve(Psi)
  U <- chol(Psi)
  d <- nrow(Psi)
  ans <- array(0, dim = c(d, d, n))
  if(!is.null(dimnames(Psi))) dimnames(ans) <- c(dimnames(Psi), list(NULL))
  ans[rep(upper.tri(Psi), n)] <- rnorm(n*d*(d-1)/2)
  ans[rep(!lower.tri(Psi, diag = FALSE) &
            !upper.tri(Psi, diag = FALSE), n)] <- sqrt(rchisq(n*d, df = nu-1:d+1))
  for(ii in 1:n) {
    tmp <- ans[,,ii] %*% U
    if(inv) tmp <- backsolve(tmp, diag(d), transpose = TRUE)
    ans[,,ii] <- crossprod(tmp)
  }
  ans
}

#################################################

#'@title Random draws from a Normal-Inverse-Wishart distribution.
#'@description Generates random draws from a Normal-Inverse-Wishart (NIW) distribution.
#'Can be used to compare prior to posterior parameter distributions.
#'@details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#'\deqn{\Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).}
#'@param n number of samples to draw.
#'@param lambda location parameter. See Details.
#'@param kappa scale parameter. See Details.
#'@param Psi scale matrix.  See Details
#'@param nu degrees of freedom.  See Details.
#'@return Returns a list with elements \code{mu} and \code{Sigma} of sizes \code{c(n,length(lambda))} and \code{c(nrow(Psi),ncol(Psi),n)}.
#'@examples
#'d <- 4 # number of dimensions
#'nu <- 7 # degrees of freedom
#'Psi <- crossprod(matrix(rnorm(d^2), d, d)) # scale
#'lambda <- rnorm(d)
#'kappa <- 2
#'n <- 1e4
#'
#'niw.sim <- rniw(n, lambda, kappa, Psi, nu)
#'
#'# diagonal elements of Sigma^{-1} are const * chi^2
#'S <- apply(niw.sim$Sigma, 3, function(M) diag(solve(M)))
#'
#'ii <- 2
#'const <- solve(Psi)[ii,ii]
#'hist(S[ii,], breaks = 100, freq = FALSE,
#'     main = parse(text = paste0("\"Histogram of \"*(Sigma^{-1})[", ii,ii,"]")),
#'     xlab = parse(text = paste0("(Sigma^{-1})[", ii,ii,"]")))
#'curve(dchisq(x/const, df = nu)/const,
#'      from = min(S[ii,]), to = max(S[ii,]), col = "red", add = TRUE)
#'
#'# elements of mu have a t-distribution
#'mu <- niw.sim$mu
#'
#'ii <- 4
#'const <- sqrt(Psi[ii,ii]/(kappa*(nu-d+1)))
#'hist(mu[,ii], breaks = 100, freq = FALSE,
#'     main = parse(text = paste0("\"Histogram of \"*mu[", ii, "]")),
#'     xlab = parse(text = paste0("mu[", ii, "]")))
#'curve(dt((x-lambda[ii])/const, df = nu-d+1)/const, add = TRUE, col = "red")
#'@seealso \code{\link{rwish}}, \code{\link{niw.mom}}, \code{\link{niw.coeffs}}.
#'@export
rniw <- function(n, lambda, kappa, Psi, nu) {
  d <- length(lambda)
  Sigma <- rwish(n, Psi, nu, inv = TRUE)
  mu <- matrix(NA, n, d)
  colnames(mu) <- names(lambda)
  for(ii in 1:n) {
    mu[ii,] <- rmvnorm(1, mean = lambda, sigma = Sigma[,,ii]/kappa)
  }
  list(mu = mu, Sigma = Sigma)
}

#################################################

#'@title Mean and variance of the Normal-Inverse-Wishart distribution.
#'@description This function computes the mean and variance of the Normal-Inverse-Wishart (NIW)
#'distribution.  Can be used to very quickly compute Bayesian point estimates for the conjugate
#'NIW prior.
#'@details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#'\deqn{\Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).}
#'Note that cov\eqn{(\mu, \Sigma) = 0}.
#'@param lambda location parameter. See Details.
#'@param kappa scale parameter. See Details.
#'@param Psi scale matrix.  See Details
#'@param nu degrees of freedom.  See Details.
#'@return Returns a list with elements \code{mu} and \code{Sigma}, each containing lists with
#'elements \code{mean} and \code{var}.  For \code{mu} these elements are of size \code{length(lambda)}
#'and \code{c(length(lambda),length(lambda))}.  For \code{Sigma} they are of size \code{dim(Psi)}
#' and \code{c(dim(Psi), dim(Psi))}, such that cov\eqn{(\Sigma_{ij}, \Sigma_{kl})=}\code{Sigma$var[i,j,k,l]}.
#'@seealso \code{\link{rniw}}, \code{\link{niw.coeffs}}, \code{\link{niw.post}}.
#'@examples
#'# NIW parameters
#'d <- 3 # number of dimensions
#'lambda <- rnorm(d)
#'kappa <- 2
#'Psi <- crossprod(matrix(rnorm(d^2), d, d))
#'nu <- 10
#'
#'# simulate data
#'niw.sim <- rniw(n = 1e4, lambda, kappa, Psi, nu)
#'
#'# check moments
#'niw.mV <- niw.mom(lambda, kappa, Psi, nu)
#'
#'# mean of mu
#'ii <- 1
#'c(true = niw.mV$mu$mean[ii], sim = mean(niw.sim$mu[,ii]))
#'
#'# variance of Sigma
#'II <- c(1,2)
#'JJ <- c(2,3)
#'c(true = niw.mV$var[II[1],II[2],JJ[1],JJ[2]],
#'  sim = cov(niw.sim$Sigma[II[1],II[2],], niw.sim$Sigma[JJ[1],JJ[2],]))
#'@export
niw.mom <- function(lambda, kappa, Psi, nu) {
  d <- length(lambda)
  b <- nu-d
  mu.mean <- lambda
  Sigma.mean <- Psi/(b-1)
  mu.var <- Sigma.mean/kappa
  Sigma.var <- Psi %o% Psi
  Sigma.var <- 2 * Sigma.var + (b-1) * (aperm(Sigma.var, c(1,3,2,4)) +
                                        aperm(Sigma.var, c(1,4,3,2)))
  Sigma.var <- Sigma.var/(b*(b-1)^2*(b-3))
  list(mu = list(mean = mu.mean, var = mu.var),
       Sigma = list(mean = Sigma.mean, var = Sigma.var))
}

#################################################

#'@title Posterior coefficients of the Normal-Inverse-Wishart distribution with its conjugate prior.
#'@description Given iid \eqn{d}-dimensional niche indicators \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)},
#'this function calculates the coefficients of the Normal-Inverse-Wishart (NIW) posterior
#'\eqn{p(\mu, \Sigma | X)} for a conjugate NIW prior.  Together with \code{\link{niw.mom}},
#'this can be used to rapidly compute the point estimates \eqn{E[\mu | X]} and \eqn{E[\Sigma | X]}.
#'@details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#'\deqn{\Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).}
#'The default value \code{kappa = 0} uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.
#'The default value \code{Psi = 0} uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'The default value \code{nu = ncol(X)+1} for \code{kappa = 0} and \code{Psi = 0} makes \eqn{E[\mu|X]=\code{colMeans(X)}} and \eqn{E[\Sigma | X]=\code{var(X)}}.
#'@param X a data matrix with observations along the rows.
#'@param lambda location parameter. See Details.
#'@param kappa scale parameter. Defaults to \code{kappa = 0}.  See Details.
#'@param Psi scale matrix. Defaults to \code{Psi = 0}.  See Details.
#'@param nu degrees of freedom. Defaults to \code{nu = ncol(X)+1}.  See Details.
#'@return Returns a list with elements \code{lambda}, \code{kappa}, \code{Psi}, \code{nu} corresponding to the coefficients of the NIW
#'posterior distribution \eqn{p(\mu, \Sigma | X)}.
#'@seealso \code{\link{rniw}}, \code{\link{niw.mom}}, \code{\link{niw.post}}.
#'@examples
#'# NIW prior coefficients
#'d <- 3
#'lambda <- rnorm(d)
#'kappa <- 5
#'Psi <- crossprod(matrix(rnorm(d^2), d, d))
#'nu <- 10
#'
#'# data
#'data(fish)
#'X <- fish[fish$species == "ARCS",2:4]
#'
#'# NIW posterior coefficients
#'post.coef <- niw.coeffs(X, lambda, kappa, Psi, nu)
#'
#'# compare
#'mu.mean <- niw.mom(post.coef$lambda, post.coef$kappa, post.coef$Psi, post.coef$nu)$mu$mean
#'mu.est <- rbind(prior = niw.mom(lambda, kappa, Psi, nu)$mu$mean,
#'                data = colMeans(X),
#'                post = mu.mean)
#'round(mu.est, 2)
#'@export
niw.coeffs <- function(X, lambda, kappa, Psi, nu) {
  d <- ncol(X)
  N <- nrow(X)
  if(missing(kappa)) kappa <- 0
  if(missing(Psi)) Psi <- 0
  if(missing(nu)) nu <- ncol(X)+1
  # sufficient statistics
  Xbar <- colMeans(X)
  S <- t(X)-Xbar
  S <- S %*% t(S)
  # posterior parameter values
  Psi2 <- Psi + S
  lambda2 <- N*Xbar
  if(kappa != 0) {
    Psi2 <- Psi2 + (N*kappa)/(N+kappa) * (Xbar-lambda) %*% t(Xbar-lambda)
    lambda2 <- lambda2 + kappa*lambda
  }
  lambda2 <- lambda2/(N+kappa)
  nu2 <- N+nu-(kappa==0)
  kappa2 <- N+kappa
  list(lambda = lambda2, kappa = kappa2, Psi = Psi2, nu = nu2)
}

#################################################

#'@title Random draws from the posterior distribution with Normal-Inverse-Wishart (NIW) prior.
#'@description Given iid \eqn{d}-dimensional niche indicators  \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)},
#'this function generates random draws from \eqn{p(\mu,\Sigma | X)} for the Normal-Inverse-Wishart
#'(NIW) prior.
#'@details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#'\deqn{\Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).}
#'The default value \code{kappa = 0} uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.
#'The default value \code{Psi = 0} uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'The default value \code{nu = ncol(X)+1} for \code{kappa = 0} and \code{Psi = 0} makes \eqn{E[\mu|X]=\code{colMeans(X)}} and \eqn{E[\Sigma | X]=\code{var(X)}}.
#'@param nsamples the number of posterior draws.
#'@param X a data matrix with observations along the rows.
#'@param lambda location parameter. See Details.
#'@param kappa scale parameter. Defaults to \code{kappa = 0}.  See Details.
#'@param Psi scale matrix. Defaults to \code{Psi = 0}.  See Details.
#'@param nu degrees of freedom. Defaults to \code{nu = ncol(X)+1}.  See Details.
#'@return Returns a list with elements \code{mu} and \code{Sigma} of sizes \code{c(nsamples, length(lambda))} and \code{c(dim(Psi), nsamples)}.
#'@seealso \code{\link{rniw}}, \code{\link{niiw.post}}.
#'@examples
#'# compare the default non-informative prior to an arbitrary informative prior
#'# for simulated data
#'
#'# simulate data
#'d <- 4
#'mu0 <- rnorm(d)
#'Sigma0 <- matrix(rnorm(d^2), d, d)
#'Sigma0 <- Sigma0 %*% t(Sigma0)
#'N <- 1e2
#'X <- rmvnorm(N, mean = mu0, sigma = Sigma0)
#'
#'# informative prior parameters
#'lambda <- rnorm(d)
#'kappa <- 20
#'Psi <- crossprod(matrix(rnorm(d^2), d, d))
#'nu <- 5
#'
#'# iid draws from informative prior pi(mu, Sigma)
#'nsamples <- 2e3
#'siw0 <- rniw(nsamples, lambda, kappa, Psi, nu)
#'
#'# iid draws from posterior p(mu, Sigma | X) with informative prior
#'siw1 <- niw.post(nsamples, X, lambda, kappa, Psi, nu)
#'
#'# iid draws from posterior p(mu, Sigma | X) with default noninformative prior
#'siw2 <- niw.post(nsamples, X)
#'
#'# compare
#'
#'# prior and posterior densities of mu
#'clrs <- c("orange", "red", "blue", "black")
#'ii <- 1
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(list(siw0, siw1, siw2), col = clrs[1:3],
#'               plot.index = ii, ylab = "Density")
#'abline(v = mu0[ii], col = clrs[4]) # true value of mu
#'legend(x = "topright",
#'       legend = c(parse(text = paste0("pi(mu[", ii, "])")),
#'                  parse(text = paste0("p(mu[", ii, "]*\" | \"*X)*\", Informative Prior\"")),
#'                  parse(text = paste0("p(mu[", ii, "]*\" | \"*X)*\", Noninformative Prior\"")),
#'                  parse(text = paste0("\"True value of \"*mu[", ii, "]"))),
#'       fill = clrs)
#'
#'# prior and posterior densities of Sigma
#'ii <- 1
#'jj <- 2
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(list(siw0, siw1, siw2), col = clrs[1:3],
#'               plot.index = c(ii,jj), ylab = "Density")
#'abline(v = Sigma0[ii,jj], col = clrs[4])
#'legend(x = "topright",
#'       legend = c(parse(text = paste0("pi(Sigma[", ii, "*", jj, "])")),
#'                  parse(text = paste0("p(Sigma[", ii, "*", jj,
#'                                      "]*\" | \"*X)*\", Informative Prior\"")),
#'                  parse(text = paste0("p(Sigma[", ii, "*", jj,
#'                                      "]*\" | \"*X)*\", Noninformative Prior\"")),
#'                  parse(text = paste0("\"True value of \"*Sigma[", ii, "*", jj, "]"))),
#'       fill = clrs)
#'
#'@export
niw.post <- function(nsamples, X, lambda, kappa, Psi, nu) {
  par <- niw.coeffs(X, lambda, kappa, Psi, nu)
  rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
}

#################################################

#'@title Random draws from the posterior distribution with Normal-Independent-Inverse-Wishart (NIIW) prior.
#'
#'@description Given iid \eqn{d}-dimensional niche indicators  \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)},
#'this function generates random draws from \eqn{p(\mu,\Sigma | X)} for the Normal-Independent-Inverse-Wishart (NIIW) prior.
#'@details The NIIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#'\deqn{\Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Omega).}
#'The default value \code{Omega = 0} uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.  In this case the NIW and NIIW priors produce identical resuls, but \code{\link{niw.post}} is faster.
#'The default value \code{Psi = 0} uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'The default value \code{nu = ncol(X)+1} for \code{Omega = 0} and \code{Psi = 0} makes \eqn{E[\mu|X]=\code{colMeans(X)}} and \eqn{E[\Sigma | X]=\code{var(X)}}.
#'Random draws are obtained by a Markov chain Monte Carlo (MCMC) algorithm; specifically,
#'a Gibbs sampler alternates between draws from \eqn{p(\mu | \Sigma, X)} and \eqn{p(\Sigma | \mu, X)}, which are Normal and Inverse-Wishart distributions respectively.
#'@param nsamples the number of posterior draws.
#'@param X a data matrix with observations along the rows.
#'@param lambda mean of mu. See Details.
#'@param Omega variance of mu. Defaults to \code{Omega = 0}.  See Details.
#'@param Psi scale matrix of Sigma. Defaults to \code{Psi = 0}.  See Details.
#'@param nu degrees of freedom of Sigma. Defaults to \code{nu = ncol(X)+1}.  See Details.
#'@param mu0 initial value of mu to start the Gibbs sampler.  See Details.
#'@param burn burn-in for the MCMC sampling algorithm.  Either an integer giving the number of initial samples to discard, or a fraction with \code{0 < burn < 1}.  Defaults to \code{burn = floor(nsamples/10)}.
#'@return Returns a list with elements \code{mu} and \code{Sigma} of sizes \code{c(nsamples, length(lambda))} and \code{c(dim(Psi), nsamples)}.
#'@seealso \code{\link{niw.post}}, \code{\link{rwish}}.
#'@examples
#'# simulate data
#'d <- 4
#'mu0 <- rnorm(d)
#'Sigma0 <- matrix(rnorm(d^2), d, d)
#'Sigma0 <- Sigma0 %*% t(Sigma0)
#'N <- 100
#'X <- rmvnorm(N, mean = mu0, sigma = Sigma0)
#'
#'# prior parameters
#'# flat prior on mu
#'lambda <- 0
#'Omega <- 0
#'# informative prior on Sigma
#'Psi <- crossprod(matrix(rnorm(d^2), d, d))
#'nu <- 5
#'
#'# sample from NIIW posterior
#'nsamples <- 2e3
#'system.time({
#'  siiw <- niiw.post(nsamples, X, lambda, Omega, Psi, nu, burn = 100)
#'})
#'
#'# sample from NIW posterior
#'kappa <- 0
#'system.time({
#'  siw <- niw.post(nsamples, X, lambda, kappa, Psi, nu)
#'})
#'
#'# check that posteriors are the same
#'
#'# p(mu | X)
#'clrs <- c("black", "red")
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(list(siiw, siw), col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
#'legend(x = "topright", legend = c("NIIW Prior", "NIW Prior"), fill = clrs)
#'
#'# p(Sigma | X)
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(list(siiw, siw), col = clrs, plot.mu = FALSE, plot.Sigma = TRUE)
#'legend(x = "topright", legend = c("NIIW Prior", "NIW Prior"), fill = clrs)
#'@export
niiw.post <- function(nsamples, X, lambda, Omega, Psi, nu, mu0 = lambda, burn) {
  # sufficient statistics
  d <- ncol(X)
  N <- nrow(X)
  Xbar <- colMeans(X)
  S <- t(X)-Xbar
  S <- S %*% t(S)
  # local variables
  mu <- rep(0, d)
  Sigma <- matrix(0, d, d)
  if(missing(burn)) burn <- .1
  if(burn < 1) burn <- floor(nsamples*burn)
  mu <- mu0
  # output
  mu.out <- matrix(NA, nsamples, d)
  Sigma.out <- array(NA, dim = c(d, d, nsamples))
  # main for loop
  for(ii in (-burn+1):nsamples) {
    # sample Sigma
    Psi2 <- N * ((mu-Xbar) %o% (mu-Xbar)) + S + Psi
    nu2 <- N+nu
    Sigma <- matrix(rwish(1, Psi2, nu2, inv = TRUE), d, d)
    # sample mu
    Sigma2 <- Sigma/N
    if(!all(Omega == 0)) {
      B <- Sigma2 %*% solve(Omega + Sigma2)
    } else B <- matrix(0, d, d)
    IB <- diag(d)-B
    lambda2 <- c(IB %*% Xbar)
    if(!all(Omega == 0)) lambda2 <- lambda2 + c(B %*% lambda)
    Omega2 <- IB %*% Sigma2
    mu <- c(rmvnorm(1, lambda2, Omega2))
    # store
    if(ii > 0) {
      mu.out[ii,] <- mu
      Sigma.out[,,ii] <- Sigma
    }
  }
  list(mu = mu.out, Sigma = Sigma.out)
}

#################################################

#'@title Plot for niche parameters
#'
#'@description For one or more species, plots some or all of the niche parameters \eqn{\mu} and \eqn{\Sigma}.
#'@param niche.par list with \code{nspecies = length(niche.par)}, each element of which is a list with parameters \code{mu} and \code{Sigma}.  See Details.
#'@param plot.mu logical.  If \code{TRUE}, plot the distribution of \eqn{\mu} for each niche indicator (e.g., stable isotope).  See Details.
#'@param plot.Sigma logical.  If \code{TRUE}, plot the distribution of \eqn{\Sigma} for each niche indicator.  See Details.
#'@param plot.index either a scalar of a numeric vector of length 2.  If \code{plot.index = i} then plot the distribution of \eqn{\mu_i}.  If \code{plot.index = c(i,j)} then plot the distribution of \eqn{\Sigma_{ij}}.
#'@param col vector of colors in which to plot each species.
#'@param ndens number of points at which to evaluate density estimates.
#'@param ylab optional label for \eqn{y}-axis.  If missing, defaults to \eqn{p(\mu_i | X)} and \eqn{p(\Sigma_{ij} | X)}.
#'@details \code{niche.par} is a list, each element of which is a distribution of niche parameters.  That is, \code{names(niche.par[[1]]) = c("mu", "Sigma")}, and if \code{niso} is the number of niche indicators (e.g., stable isotopes), then \code{dim(niche.par[[1]]$mu) = c(nsamples, niso)} and \code{dim(niche.par[[1]]$Sigma) = c(niso, niso, nsamples)}.
#'@seealso \code{\link{niw.post}}, \code{\link{niiw.post}} for niche parameter output, \code{density} in the \code{R} \code{base} package for density estimation from sample data.
#'@return Returns a plot of the distribution of some or all niche parameters.
#'@examples
#'# fish data
#'data(fish)
#'
#'# generate parameter draws from the "default" posteriors of each fish
#'nsamples <- 1e3
#'system.time({
#'  fish.par <- tapply(1:nrow(fish), fish$species,
#'                     function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
#'})
#'
#'# various parameter plots
#'clrs <- c("black", "red", "blue", "orange") # colors for each species
#'
#'# mu1, mu2, and Sigma12
#'par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
#'niche.par.plot(fish.par, col = clrs, plot.index = 1)
#'niche.par.plot(fish.par, col = clrs, plot.index = 2)
#'niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
#'legend("topright", legend = names(fish.par), fill = clrs)
#'
#'# all mu
#'niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
#'legend("topright", legend = names(fish.par), fill = clrs)
#'
#'# all mu and Sigma
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
#'legend("topright", legend = names(fish.par), fill = clrs)
#'@export
niche.par.plot <- function(niche.par, plot.mu = TRUE, plot.Sigma = TRUE, plot.index,
                           col, ndens = 512, ylab) {
  niso <- ncol(niche.par[[1]]$mu)
  nsmp <- length(niche.par)
  # determine the number of rows and columns in plot
  if(!missing(plot.index)) {
    if(length(plot.index) == 1) {
      plot.mu <- TRUE
      plot.Sigma <- FALSE
    } else if(length(plot.index) == 2) {
      plot.mu <- FALSE
      plot.Sigma <- TRUE
    } else {
      stop("Incorrect specification of plot.index.  Must be a numeric vector of length 1 (plot.mu) or 2 (plot.Sigma).")
    }
    nc <- 1
    nr <- 1
  } else {
    nc <- niso
    nr <- 0
    if(plot.mu) nr <- nr + 1
    if(plot.Sigma) nr <- nr + niso
  }
  # determine ylab
  if(missing(ylab)) {
    if(plot.mu) {
      ylab.mu <- sapply(1:niso, function(ii)
                        parse(text = paste0("p(mu[", ii, "]*\" | \"*X)")))
    }
    if(plot.Sigma) {
      ylab.Sigma <- apply(as.matrix(expand.grid(1:niso, 1:niso))[,2:1], 1, function(II) {
        paste0(II[1], "*", II[2])
      })
      ylab.Sigma <- sapply(ylab.Sigma, function(ii) {
        parse(text = paste0("p(Sigma[", ii, "]*\" | \"*X)"))
      })
    }
  } else {
    ylab <- rep(ylab, len = niso*(niso+1))
    ylab.mu <- ylab[1:niso]
    ylab.Sigma <- ylab[-(1:niso)]
  }
  # plot
  densx <- matrix(NA, ndens, nsmp)
  densy <- densx
  if((nr > 1) || (nc > 1)) par(mfrow = c(nr, nc))
  # mu
  if(plot.mu) {
    for(ii in 1:niso) {
      if(missing(plot.index) || plot.index == ii) {
        for(kk in 1:nsmp) {
          dens <- density(niche.par[[kk]]$mu[,ii], n = ndens)
          densx[,kk] <- dens$x
          densy[,kk] <- dens$y
        }
        plot(densx, densy, type = "n",
             xlab = parse(text = paste0("mu[", ii, "]")),
             ylab = ylab.mu[ii])
        for(kk in 1:nsmp) lines(densx[,kk], densy[,kk], col = col[kk])
      }
    }
  }
  if(plot.Sigma) {
    for(ii in 1:niso) {
      for(jj in 1:niso) {
        if(missing(plot.index) || all(plot.index == c(ii,jj))) {
          for(kk in 1:nsmp) {
            dens <- density(niche.par[[kk]]$Sigma[ii,jj,], n = ndens)
            densx[,kk] <- dens$x
            densy[,kk] <- dens$y
          }
          plot(densx, densy, type = "n",
               xlab = parse(text = paste0("Sigma[", ii, "*", jj, "]")),
               ylab = ylab.Sigma[(ii-1)*niso + jj])
          for(kk in 1:nsmp) lines(densx[,kk], densy[,kk], col = col[kk])
        }
      }
    }
  }
}

#################################################

#'@title Plot for 2-d projection of niche regions
#'
#'@description For one or more species, creates a series of plots: i) the raw niche indicators (e.g., stable isotope) data,
#'ii) their density estimates, and iii) 2-dimensional projections of probabilistic niche regions based on \eqn{n}-dimensionsional data.
#'@details A set of plots is created for each pairwise combination of niche indicators.
#'Below the diagonal are scatterplots for each species, above the diagonal are ellipses corresponding to 2-d projections of the probabilistic niche regions.  The diagonal displays density estimates for each indicator, and optionally the raw 1-d data.
#'See Swanson et al. (2014) for detailed description of the probabilistic niche region.
#'@param niche.par a list of length \code{nspecies}, each element of which in turn is a list with elements \code{mu} and \code{Sigma}.  Each of these will correspond to an ellipse being drawn for that species in the corresponding 2-d plane. See Example.
#'@param niche.data a list of length \code{nspecies}, each element of which is a matrix with observations along the rows and niche indicators (e.g., stable isotopes) along the columns.
#'@param alpha size of the niche region to plot. Defaults to 0.95.
#'@param species.names names of the species. Defaults to \code{names(niche.par)}.
#'@param iso.names names of the niche indicators (or isotopes) to plot. Defaults to \code{colnames(niche.par)}.
#'@param col vector of colours in which each species will be drawn.
#'@param ndens number of points at which to evaluate kernel density estimates.
#'@param pfrac fraction of the plot on which to display 1-dimensional raw niche indicator data. \code{pfrac = 0} means don't display the raw data in 1-d.
#'@param xlab title of plot, located on the bottom.  Defaults to no title.
#'@return Returns a series of plots displaying niche indicator data and their probabilistic niche projections.
#'@references Heidi K. Swanson, Martin Lysy, Ashley D. Stasko, Michael Power, Jim D. Johnson, and James D. Reist (2014).  ``What Would Hutchinson Think?  A Probabilistic Quantification of Multidimensional Ecological Niches and Niche Overlap''.  \emph{Ecology: Statistical Reports} (accepted).
#'@seealso \code{\link{overlap.plot}}, \code{\link{niw.post}}, \code{\link{niiw.post}}.
#'@examples
#'data(fish) # 4 fish, 3 isotopes
#'
#'# generate 10 parameter draws from the posteriors of each fish with default prior
#'nsamples <- 10
#'fish.par <- tapply(1:nrow(fish), fish$species,
#'                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
#'
#'# format data for plotting function
#'fish.data <- tapply(1:nrow(fish), fish$species, function(ii) X = fish[ii,2:4])
#'
#'clrs <- c("black", "red", "blue", "orange") # colors for each species
#'niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .1,
#'           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
#'           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))
#'@export
niche.plot <- function(niche.par, niche.data, alpha = .95,
                       species.names, iso.names,
                       col, ndens = 512, pfrac = 0, xlab) {
  niso <- ncol(niche.par[[1]]$mu)
  nspec <- length(niche.par)
  npts <- 100 # number of points for each ellipse
  nell <- sapply(niche.par, function(x) nrow(x$mu)) # number of ellipses per species
  if(missing(species.names)) species.names <- names(niche.par)
  if(missing(iso.names)) iso.names <- colnames(niche.par[[1]]$mu)
  # create all the ellipses to get the plot limits right.
  ell <- vector("list", nspec)
  names(ell) <- names(species.names)
  D <- combn(niso, 2)
  for(ii in 1:nspec) {
    ell.tmp <- array(NA, c(nell[ii], ncol(D), npts+1, 2))
    for(jj in 1:nell[ii]) {
      for(kk in 1:ncol(D)) {
        ell.coord <- ellipse(niche.par[[ii]]$mu[jj, D[,kk]],
                             V = niche.par[[ii]]$Sigma[D[,kk], D[,kk], jj],
                             alpha = alpha, n = npts)
        ell.tmp[jj,kk,,] <- ell.coord
      }
    }
    ell[[ii]] <- ell.tmp
  }
  # plot limits.
  lims <- array(sapply(niche.data, function(x) apply(x, 2, range)),
                dim = c(2, niso, nspec))
  lims <- apply(lims, 2, range)
  # plots
  par(mfcol = c(niso,niso), mar = rep(.5, 4), oma = rep(4,4))
  for(ci in 1:niso) {
    for(ri in 1:niso) {
      # initialize plot
      plot.new()
      plot.window(lims[,ci], lims[,ri])
      if (ci == ri) {
        # diagonals: density plots
        xdens <- matrix(NA, ndens, nspec)
        ydens <- xdens
        for(ii in 1:nspec) {
          den <- density(niche.data[[ii]][,ci], n = ndens)
          xdens[,ii] <- den$x
          ydens[,ii] <- den$y
        }
        for(ii in 1:nspec) {
          ly <- par("usr")[1:2]
          ly[2] <- ly[1] + pfrac*(ly[2]-ly[1])
          ly[3] <- (ly[2]-ly[1])/nspec
          segments(x0 = niche.data[[ii]][,ci],
                   y0 = ly[1]+(ii-1)*ly[3], y1 = ly[1]+ii*ly[3], col = col[ii])
          ly <- ly[2] + ydens[,ii]/max(ydens)*(lims[2,ci]-ly[2])
          lines(xdens[,ii], ly, col = col[ii])
        }
      }
      if (ri > ci) {
        # lower triangle: point plots
        for(ii in 1:nspec) {
          points(niche.data[[ii]][,c(ci,ri)], col = col[ii], pch = 16)
        }
      }
      if (ri < ci) {
        # upper triangle: ellipses
        for(ii in 1:nspec) {
          for(jj in 1:nell[ii]) {
            lines(ell[[ii]][jj,which(D[1,] == ri & D[2,] == ci),,2:1], col = col[ii])
          }
        }
      }
      box()
      if(ci == niso) axis(side = 4) else axis(side = 4, labels = FALSE)
      if(ri == niso) axis(side = 1) else axis(side = 1, labels = FALSE)
      if(ri == 1) mtext(text = iso.names[ci], side = 3, line = 1)
      if(ci == 1) mtext(text = iso.names[ri], side = 2, line = 1)
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, outer = TRUE, line = 2.2, cex = .9)
  }
  legend(x = "topleft", legend = species.names, fill = col, bty = "n", cex = 1.25)
}

#################################################

#'@title Monte Carlo calculation of niche region overlap metrics
#'
#'@description Calculates the distribution of a niche region overlap metric for each pairwise species combination
#'and user-specified niche region sizes. 
#'@details The overlap metric is the probability that a randomly drawn individual from
#'species \eqn{A} will be found within the niche region of species \eqn{B} (for a given niche region 
#'size, e.g., \code{alpha = .95}).  It is a single number which is a function of the
#'parameters for each species, \eqn{\Theta_A = (\mu_A, \Sigma_A)} and \eqn{\Theta_B = (\mu_B, \Sigma_B)}.  This number is difficult to calculate directly, but easy to approximate stochastically by generating \code{nprob}
#'draws from the distribution of species \eqn{A} and counting the fraction of them which fall in the niche region of species \eqn{B}.
#'Typically the true values of \eqn{\Theta_A} and \eqn{\Theta_B} are unknown and must be estimated from the data.
#'Thus, the overlap metric is calculated for \code{nreps} combinations of samples from
#'\eqn{p(\Theta_A | X)} and \eqn{p(\Theta_B | X)} which are supplied in \code{niche.par}.
#'See Swanson et al. (2014) for a detailed description of niche overlap and its calculation.
#'@param niche.par a list with \code{nspecies = length(niche.par)}, each element of which in turn is a list with elements \code{mu} and \code{Sigma}.  See Details.
#'@param nreps the number of overlap metric calculations for each species.  Defaults to
#'the smallest number of parameter samples supplied by \code{niche.par}.  See Details.
#'@param nprob the number of normal draws for each Monte Carlo overlap metric calculation.  See Details.
#'@param alpha scalar or vector of niche region sizes for calculating the niche overlap metric. Defaults to 0.95.
#'@param species.names names of the species. Defaults to \code{names(niche.par)}.
#'@param norm.redraw logical. If \code{FALSE}, the same \code{nprob*nspecies} iid \eqn{N(0,1)}
#'draws are used for each calculation of the overlap metric. This increases the Monte Carlo
#'error, but the procedure is about 1.5x faster.  Defaults to \code{TRUE}.
#'@return Returns an array of size \code{c(nspecies, nspecies, nreps, nlevels)},
#'where \code{nlevels} is the number of alpha levels at which to calculate the overlap metric.  For each of the last two dimensions of the output array, the first two dimensions form an \code{nspecies} by \code{nspecies} matrix giving
#'each pairwise calculation of overlap metric between two species for given \eqn{\Theta_A}, \eqn{\Theta_B}, and \code{alpha}.
#'In each of these matrices, Species \eqn{A} is along the rows of this matrix and Species \eqn{B} is along the columns.
#'@references Heidi K. Swanson, Martin Lysy, Ashley D. Stasko, Michael Power, Jim D. Johnson, and James D. Reist (2014).  ``What Would Hutchinson Think?  A Probabilistic Quantification of Multidimensional Ecological Niches and Niche Overlap''.  \emph{Ecology: Statistical Reports} (accepted).
#'@seealso \code{\link{overlap.plot}}, \code{\link{niw.post}}, \code{\link{niiw.post}}.
#'@examples
#'# fish data
#'data(fish)
#'
#'# generate parameter draws from the "default" posteriors of each fish
#'nsamples <- 500
#'system.time({
#'  fish.par <- tapply(1:nrow(fish), fish$species,
#'                     function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
#'})
#'
#'# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.
#'system.time({
#'  over <- overlap(fish.par, nreps = nsamples, nprob = nsamples, alpha = c(.95, .99))
#'})
#'
#'# posterior expectations of overlap metrics
#'over.mean <- apply(over*100, c(1:2, 4), mean)
#'round(over.mean)
#'
#'# posterior 95% credible intervals of overlap metrics
#'over.cred <- apply(over*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
#'round(over.cred[,,,1]) # display alpha = .95 niche region
#'@export
overlap <- function(niche.par, nreps, nprob, alpha = 0.95,
                    species.names, norm.redraw = TRUE) {
  niso <- ncol(niche.par[[1]]$mu) # number of isotopes
  nspec <- length(niche.par) # number of species
  nlevels <- length(alpha) # number of levels
  if(missing(species.names)) species.names <- names(niche.par)
  # temporary variables
  mu <- matrix(NA, niso, nspec)
  x <- array(NA, dim = c(niso, nprob, nspec))
  C <- array(NA, dim = c(niso, niso, nspec))
  Zsq <- array(NA, dim = c(nprob, nspec-1, nlevels))
  # constants
  qlevels <-  qchisq(alpha, df = niso)
  qlevels <- array(rep(qlevels, each = nprob*(nspec-1)),
                   dim = c(nprob, nspec-1, nlevels))
  if(!norm.redraw) z <- matrix(rnorm(niso*nprob), niso, nprob)
  # output
  over <- array(NA, dim = c(nspec, nspec, nreps, nlevels))
  dimnames(over) <- list(species.names, species.names, NULL,
                         paste0(alpha*100, "%"))
  names(dimnames(over)) <- c("Species A", "Species B", "", "alpha")
  # subsample the parameters posteriors of each niche nreps times
  ind <- sapply(niche.par, function(mv) {
    sample(nrow(mv$mu), nreps, replace = TRUE)
  })
  # calculate each overlap
  for(jj in 1:nreps) {
    # get means, t(chol(variance)), and draws for each species
    for(ii in 1:nspec) {
      mu[,ii] <- niche.par[[ii]]$mu[ind[jj,ii],]
      C[,,ii] <- t(chol(niche.par[[ii]]$Sigma[,,ind[jj,ii]]))
      if(norm.redraw) z <- matrix(rnorm(niso*nprob), niso, nprob)
      x[,,ii] <- C[,,ii] %*% z + mu[,ii]
    }
    # check whether the draw of each species is within the ellipse of every other
    for(ii in 1:nspec) {
      Z <- backsolve(C[,,ii], matrix(x[,,-ii], nrow = niso)-mu[,ii], upper.tri = FALSE)
      Zsq[] <- colSums(Z^2)
      over[-ii,ii,jj,] <- colMeans(Zsq < qlevels)
    }
  }
  if(dim(over)[4] == 1) over <- over[,,,1]
  over
}

#################################################

#'@title Plot the overlap metric
#'
#'@description Plots the posterior distribution of the niche region overlap metric calculated for
#'each pairwise combination of species.
#'@details This function uses the overlap metric information in \code{over.stat} calculated by \code{\link{overlap}}
#'to create 2-dimensional plots of interspecific niche region overlap.
#'@param over.stat an array with \code{dim(over.stat) = c(nspecies, nspecies, nreps)} containing \code{nreps} calculations
# of the overlap metric for each pair of species. See Details.
#'@param nbreaks number of breaks in the histogram. Defaults to 50.
#'@param equal.axis logical. If \code{TRUE}, all histograms in a given column of the output (corresponding to different Species \eqn{A} for the same Species \eqn{B}) are plotted on the same range.
#'@param species.names a vector of species names. Defaults to \code{dimnames(over.stat)[[1]]}.
#'@param col a vector of the colours in which each species will be drawn.
#'@param mean.cred logical. If \code{TRUE}, vertical lines for mean and 95\% credible intervals will be
#'included in the historgram of each overlap metric.
#'@param mean.cred.col colour of the mean and credible interval lines in the histogram.
#'@param xlab optional plot title, located on the bottom.  Default is no title.
#'@seealso \code{\link{overlap}}, \code{\link{niw.post}}, \code{\link{niiw.post}}.
#'@return Returns a series of histograms illustrating the probability of pairwise species overlap.
#'@examples
#'# fish data
#'data(fish)
#'
#'# parameter draws from the "default" posteriors of each fish
#'nsamples <- 500
#'system.time({
#'  fish.par <- tapply(1:nrow(fish), fish$species,
#'                     function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
#'})
#'
#'# overlap calculation
#'system.time({
#'  over <- overlap(fish.par, nreps = nsamples, nprob = nsamples, alpha = c(.95, .99))
#'})
#'
#'# overlap plot
#'clrs <- c("black", "red", "blue", "orange") # color for each species
#'ii <- 1 # which niche region alpha level to use
#'overlap.plot(over[,,,ii], col = clrs, mean.cred.col = "turquoise",
#'             xlab = paste0("Overlap Probability (%) -- Niche Region Size: ",
#'                           dimnames(over)[[4]][ii]))
#'@export
overlap.plot <- function(over.stat, nbreaks = 50, equal.axis = FALSE, species.names, col,
                         mean.cred = TRUE, mean.cred.col = "green", xlab) {
  if(length(dim(over.stat)) != 3 || dim(over.stat)[1] != dim(over.stat)[2])
    stop("Incorrect specification of over.stat.")
  nspec <- dim(over.stat)[1]
  if(missing(species.names)) species.names <- dimnames(over.stat)[[1]]
  # histograms
  over.hist <- apply(over.stat, 1:2,
                     function(x) {
                       if(any(is.na(x))) return(NULL)
                       else {
                         tmp <- hist(x*100, breaks = nbreaks, plot = FALSE)
                         tmp$density <- tmp$density/max(tmp$density)
                         tmp$counts <- tmp$counts/max(tmp$counts)
                       }
                       tmp
                     })
  # x-axis limits
  xlim <- lapply(over.hist,
                 function(h) {
                   if(is.null(h)) tmp <- matrix(NA, 2, 2)
                   else {
                     tmp <- cbind(range(h$breaks), range(h$density))
                   }
                   tmp
                 })
  xlim <- matrix(xlim, nspec, nspec)
  if(equal.axis) {
    for(ii in 1:nspec) {
      xlim[,ii] <- rep(list(cbind(range(sapply(xlim[,ii],
                                               function(ll) ll[,1]), na.rm = TRUE),
                                  range(sapply(xlim[,ii],
                                               function(ll) ll[,2]), na.rm = TRUE))),
                       nspec)
    }
  }
  # mean and credible intervals
  if(mean.cred) {
    over.mean <- apply(over.stat, 1:2, mean)*100
    over.cred <- apply(over.stat, 1:2, quantile, prob = c(.025, .975), na.rm = TRUE)*100
    over.cred <- array(over.cred, dim = c(2, nspec, nspec))
  }
  # plot
  par(mfcol = c(nspec,nspec), mar = c(1.5,rep(.5, 3)), oma = rep(4,4))
  for(ci in 1:nspec) {
    for(ri in 1:nspec) {
      # initialize plot
      plot.new()
      if (ri != ci) {
        # off diagonals: overlap histograms
        plot.window(xlim[ri,ci][[1]][,1], xlim[ri,ci][[1]][,2])
        if(equal.axis) {
          # recalculate histograms with new limits
          tmp <- hist(over.stat[ri,ci,]*100,
                      breaks = seq(xlim[ri,ci][[1]][1,1],
                        xlim[ri,ci][[1]][2,1], len = nbreaks+1),
                      plot = FALSE)
          tmp$density <- tmp$density/max(tmp$density)
          over.hist[[ri,ci]] <- tmp
        }
        plot(over.hist[[ri,ci]], add = TRUE, freq = FALSE, col = col[ci],
             border = "white")
        if(mean.cred) {
          # mean and 95% credible intervals
          abline(v = c(over.mean[ri,ci], over.cred[,ri,ci]),
                 col = mean.cred.col, lty = c(1,2,2), lwd = 2)
        }
      } else {
        # diagonals: empty plots
        plot.window(xlim = c(0,1), ylim = c(0,1))
      }
      if(ri == 1 && ci == 1) {
        text(x = .5, y = .5,
             labels = expression("Niche Overlap: "*bgroup("(",
                                                          atop("Row", "Column"), ")")),
             adj = c(.5, .5), cex = 1)
      }
      box()
      if(ci != ri) axis(side = 1)
      if(ci > 1) axis(side = 2, labels = FALSE)
      if(ci < nspec) axis(side = 4, labels = FALSE)
      if(ri == 1) mtext(text = species.names[ci], side = 3, line = 1, col = col[ci])
      if(ci == 1) mtext(text = species.names[ri], side = 2, line = 1)
      if(mean.cred && ri == nspec && ci == nspec) {
        legend(x = "center", legend = c("Mean", "95% Credible Interval"),
               lty = c(1, 2), lwd = 2, col = mean.cred.col)
      }
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, line = 1.5, cex = 1, outer = TRUE)
  }
}

#################################################
