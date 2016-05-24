#' Simulation of the (Restricted) Likelihood Ratio Statistic
#'
#' These functions simulate values from the (exact) finite sample distribution
#' of the (restricted) likelihood ratio statistic for testing the presence of
#' the variance component (and restrictions of the fixed effects) in a simple
#' linear mixed model with known correlation structure of the random effect and
#' i.i.d. errors. They are usually called by \code{exactLRT} or
#' \code{exactRLRT}.
#'
#' The model under the alternative must be a linear mixed model
#' \eqn{y=X\beta+Zb+\varepsilon}{y=X*beta+Z*b+epsilon} with a single random
#' effect \eqn{b} with known correlation structure \eqn{Sigma} and i.i.d errors.
#' The simulated distribution of the likelihood ratio statistic was derived by
#' Crainiceanu & Ruppert (2004). The simulation algorithm uses a gridsearch over
#' a log-regular grid of values of
#' \eqn{\lambda=\frac{Var(b)}{Var(\varepsilon)}}{lambda=Var(b)/Var(epsilon)} to
#' maximize the likelihood under the alternative for \code{nsim} realizations of
#' \eqn{y} drawn under the null hypothesis. \code{log.grid.hi} and
#' \code{log.grid.lo} are the lower and upper limits of this grid on the log
#' scale. \code{gridlength} is the number of points on the grid.\ These are just
#' wrapper functions for the underlying C code.
#'
#' @aliases RLRTSim
#' @param X The fixed effects design matrix of the model under the alternative
#' @param Z The random effects design matrix of the model under the alternative
#' @param q The number of parameters restrictions on the fixed effects (see
#'   Details)
#' @param sqrt.Sigma The upper triangular cholesky factor of the correlation
#'   matrix of the random effect
#' @param seed Specify a seed for \code{set.seed}
#' @param nsim Number of values to simulate
#' @param log.grid.hi Lower value of the grid on the log scale. See
#'   \bold{Details}
#' @param log.grid.lo Lower value of the grid on the log scale. See
#'   \bold{Details}
#' @param gridlength Length of the grid for the grid search over lambda. See
#'   \bold{Details}
#' @param parallel The type of parallel operation to be used (if any). If
#'   missing, the default is "no parallelization").
#' @param ncpus integer: number of processes to be used in parallel operation:
#'   typically one would chose this to the number of available CPUs. Defaults to
#'   1, i.e., no parallelization.
#' @param cl An optional parallel or snow cluster for use if parallel = "snow".
#'   If not supplied, a cluster on the local machine is created for the duration
#'   of the call.
#' @return A vector containig the the simulated values of the (R)LRT under the
#'   null, with attribute 'lambda' giving \eqn{\arg\min(f(\lambda))} (see
#'   Crainiceanu, Ruppert (2004)) for the simulations.
#' @author Fabian Scheipl; parallelization code adapted from \code{boot} package
#' @seealso \code{\link{exactLRT}}, \code{\link{exactRLRT}} for tests
#' @references Crainiceanu, C. and Ruppert, D. (2004) Likelihood ratio tests in
#'   linear mixed models with one variance component, \emph{Journal of the Royal
#'   Statistical Society: Series B},\bold{66},165--185.
#'
#'   Scheipl, F. (2007) Testing for nonparametric terms and random effects in
#'   structured additive regression.  Diploma thesis.\
#'   \url{http://www.statistik.lmu.de/~scheipl/downloads/DIPLOM.zip}.
#'
#'   Scheipl, F., Greven, S. and Kuechenhoff, H (2008) Size and power of tests
#'   for a zero random effect variance or polynomial regression in additive and
#'   linear mixed models, \emph{Computational Statistics & Data Analysis},
#'   \bold{52}(7):3283-3299
#' @keywords datagen distribution
#' @examples
#'
#' library(lme4)
#' g <- rep(1:10, e = 10)
#' x <- rnorm(100)
#' y <- 0.1 * x + rnorm(100)
#' m <- lmer(y ~ x + (1|g), REML=FALSE)
#' m0 <- lm(y ~ 1)
#'
#' (obs.LRT <- 2*(logLik(m)-logLik(m0)))
#' X <- getME(m,"X")
#' Z <- t(as.matrix(getME(m,"Zt")))
#' sim.LRT <- LRTSim(X, Z, 1, diag(10))
#' (pval <- mean(sim.LRT > obs.LRT))
#'
#' @export LRTSim
#' @useDynLib RLRsim
LRTSim <- function(X,Z,q, sqrt.Sigma, seed=NA, nsim=10000,
  log.grid.hi=8, log.grid.lo=-10, gridlength=200,
  parallel = c("no", "multicore", "snow"),
  ncpus = 1L, cl = NULL){
  
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore")
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow")
      have_snow <- TRUE
    if (!have_mc && !have_snow)
      ncpus <- 1L
  }
  
  K <- NCOL(Z)     # no. of random effects
  n <- NROW(X)     # no. of obs
  p <- NCOL(X)     # no of fixed effects
  
  #compute eigenvalues
  mu <- (svd(sqrt.Sigma %*% t(qr.resid(qr(X), Z)), nu = 0, nv = 0)$d)^2
  xi <- (svd(sqrt.Sigma %*% t(Z), nu = 0, nv = 0)$d)^2
  
  #norm eigenvalues
  mu <- mu / max(mu,xi)
  xi <- xi / max(mu,xi)
  
  lambda.grid <-c(0, exp(seq(log.grid.lo, log.grid.hi, length = gridlength - 1)))
  
  if (!is.na(seed))
    set.seed(seed)
  
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    nsim. <- as.integer(ceiling(nsim/ncpus))
    if (have_mc) {
      tmp <- parallel::mclapply(seq_len(ncpus), function(i){
        RLRsimCpp(p = as.integer(p), k = as.integer(K),
          n = as.integer(n), nsim = as.integer(nsim.),
          g = as.integer(gridlength),
          q = as.integer(q), mu = as.double(mu),
          lambda = as.double(lambda.grid),
          lambda0 = as.double(0), xi = as.double(xi),
          REML = as.logical(FALSE))
      }, mc.cores = ncpus)
      do.call(mapply, c(tmp, FUN=c))
    } else {
      if (have_snow) {
        if (is.null(cl)) {
          cl <- parallel::makePSOCKcluster(rep("localhost",
            ncpus))
          if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
            parallel::clusterSetRNGStream(cl)
          }
          tmp <- parallel::parLapply(cl, seq_len(ncpus), function(i){
            RLRsimCpp(p = as.integer(p), k = as.integer(K),
              n = as.integer(n), nsim = as.integer(nsim.),
              g = as.integer(gridlength),
              q = as.integer(q), mu = as.double(mu),
              lambda = as.double(lambda.grid),
              lambda0 = as.double(0), xi = as.double(xi),
              REML = as.logical(FALSE))
          })
          parallel::stopCluster(cl)
          do.call(mapply, c(tmp, FUN=c))
        } else {
          tmp <- parallel::parLapply(cl, seq_len(ncpus), function(i){
            RLRsimCpp(p = as.integer(p), k = as.integer(K),
              n = as.integer(n), nsim = as.integer(nsim.),
              g = as.integer(gridlength),
              q = as.integer(q), mu = as.double(mu),
              lambda = as.double(lambda.grid),
              lambda0 = as.double(0), xi = as.double(xi),
              REML = as.logical(FALSE))
          })
          do.call(mapply, c(tmp, FUN=c))
        }
      }
    }
  } else {
    RLRsimCpp(p = as.integer(p), k = as.integer(K),
      n = as.integer(n), nsim = as.integer(nsim),
      g = as.integer(gridlength),
      q = as.integer(q), mu = as.double(mu),
      lambda = as.double(lambda.grid),
      lambda0 = as.double(0), xi = as.double(xi),
      REML = as.logical(FALSE))
  }
  lambda <- lambda.grid[res$lambdaind]
  ret <- res$res
  attr(ret, "lambda") <- lambda.grid[res$lambdaind]
  return(ret)
}

