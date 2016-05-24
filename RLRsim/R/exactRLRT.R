#' Restricted Likelihood Ratio Tests for additive and linear mixed models
#' 
#' This function provides an (exact) restricted likelihood ratio test based on
#' simulated values from the finite sample distribution for testing whether the
#' variance of a random effect is 0 in a linear mixed model with known
#' correlation structure of the tested random effect and i.i.d. errors.
#' 
#' Testing in models with only a single variance component require only the
#' first argument \code{m}. For testing in models with multiple variance
#' components, the fitted model \code{m} must contain \bold{only} the random
#' effect set to zero under the null hypothesis, while \code{mA} and \code{m0}
#' are the models under the alternative and the null, respectively. For models
#' with a single variance component, the simulated distribution is exact if the
#' number of parameters (fixed and random) is smaller than the number of
#' observations. Extensive simulation studies (see second reference below)
#' confirm that the application of the test to models with multiple variance
#' components is safe and the simulated distribution is correct as long as the
#' number of parameters (fixed and random) is smaller than the number of
#' observations and the nuisance variance components are not superfluous or
#' very small. We use the finite sample distribution of the restricted
#' likelihood ratio test statistic as derived by Crainiceanu & Ruppert (2004).
#' 
#' @param m The fitted model under the alternative or, for testing in models
#' with multiple variance components, the reduced model containing only the
#' random effect to be tested (see Details), an \code{lme}, \code{lmerMod} or
#' \code{spm} object
#' @param mA The full model under the alternative for testing in models with
#' multiple variance components
#' @param m0 The model under the null for testing in models with multiple
#' variance components
#' @param seed input for \code{set.seed}
#' @param nsim Number of values to simulate
#' @param log.grid.hi Lower value of the grid on the log scale. See
#' \code{\link{exactRLRT}}.
#' @param log.grid.lo Lower value of the grid on the log scale. See
#' \code{\link{exactRLRT}}.
#' @param gridlength Length of the grid. See \code{\link{exactLRT}}.
#' @param parallel The type of parallel operation to be used (if any). If
#' missing, the default is "no parallelization").
#' @param ncpus integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs. Defaults to
#' 1, i.e., no parallelization.
#' @param cl An optional parallel or snow cluster for use if parallel = "snow".
#' If not supplied, a cluster on the local machine is created for the duration
#' of the call.
#' @return A list of class \code{htest} containing the following components:
#' @return A list of class \code{htest} containing the following components:
#' \itemize{
#' \item \code{statistic} the observed likelihood ratio
#' \item \code{p} p-value for the observed test statistic
#' \item \code{method} a character string indicating what type of test was
#' performed and how many values were simulated to determine the critical value
#' \item \code{sample} the samples from the null distribution returned by
#' \code{\link{RLRTSim}}
#' }
#' @author Fabian Scheipl, bug fixes by Andrzej Galecki, updates for
#' \pkg{lme4}-compatibility by Ben Bolker
#' @seealso \code{\link{RLRTSim}} for the underlying simulation algorithm;
#' \code{\link{exactLRT}} for likelihood based tests
#' @references Crainiceanu, C. and Ruppert, D. (2004) Likelihood ratio tests in
#' linear mixed models with one variance component, \emph{Journal of the Royal
#' Statistical Society: Series B},\bold{66},165--185.
#' 
#' Greven, S., Crainiceanu, C., Kuechenhoff, H., and Peters, A. (2008)
#' Restricted Likelihood Ratio Testing for Zero Variance Components in Linear
#' Mixed Models, \emph{Journal of Computational and Graphical Statistics},
#' \bold{17} (4): 870--891.
#' 
#' Scheipl, F., Greven, S. and Kuechenhoff, H. (2008) Size and power of tests
#' for a zero random effect variance or polynomial regression in additive and
#' linear mixed models.  \emph{Computational Statistics & Data Analysis},
#' \bold{52}(7):3283--3299.
#' @keywords htest
#' @examples
#' 
#' library(lme4)
#' data(sleepstudy)
#' mA <- lmer(Reaction ~ I(Days-4.5) + (1|Subject) + (0 + I(Days-4.5)|Subject), 
#'   data = sleepstudy)
#' m0 <- update(mA, . ~ . - (0 + I(Days-4.5)|Subject))
#' m.slope  <- update(mA, . ~ . - (1|Subject))
#' #test for subject specific slopes:
#' exactRLRT(m.slope, mA, m0)
#' 
#' library(mgcv)
#' data(trees)
#' #test quadratic trend vs. smooth alternative
#' m.q<-gamm(I(log(Volume)) ~ Height + s(Girth, m = 3), data = trees, 
#'   method = "REML")$lme
#' exactRLRT(m.q)
#' #test linear trend vs. smooth alternative
#' m.l<-gamm(I(log(Volume)) ~ Height + s(Girth, m = 2), data = trees, 
#'   method = "REML")$lme
#' exactRLRT(m.l)
#' 
#' @export exactRLRT
#' @importFrom stats anova cov2cor logLik quantile 
'exactRLRT' <- function(m, mA = NULL, m0 = NULL, seed = NA, 
  nsim = 10000, log.grid.hi = 8, log.grid.lo = -10, gridlength = 200,
  parallel = c("no", "multicore", "snow"), 
  ncpus = 1L, cl = NULL) {
  if (class(m) == "spm") {
    m <- m$fit
    class(m) <- "lme"
  }
  if (class(m) %in% c("amer", "mer")) 
    stop("Models fit with package <amer> or versions of <lme4> below 1.0 are no longer supported.")
  if (!(c.m <- (class(m))) %in% c("lme", "lmerMod", "merModLmerTest")) 
    stop("Invalid <m> specified. \n")
  if ("REML" != switch(c.m, 
    lme = m$method, 
    lmerMod = ifelse(lme4::isREML(m), "REML", "ML"))){
    message("Using restricted likelihood evaluated at ML estimators.")
    message("Refit with method=\"REML\" for exact results.")	
  }
  
  d <- switch(c.m, lme = extract.lmeDesign(m), 
    lmerMod = extract.lmerModDesign(m))
  X <- d$X
  qrX <- qr(X)
  Z <- d$Z
  y <- d$y
  Vr <- d$Vr
  if(all(Vr == 0)){
    # this only happens if the estimate of the tested variance component is 0. 
    # since we still want chol(cov2cor(Vr)) to work, this does the trick.
    diag(Vr) <- 1
  }
  K <- ncol(Z)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(mA) && is.null(m0)) {
    if(length(d$lambda) != 1 || d$k != 1) 
      stop("multiple random effects in model - 
                 exactRLRT needs <m> with only a single random effect.")
    #2*restricted ProfileLogLik under H0: lambda=0
    res <- qr.resid(qrX, y)
    R <- qr.R(qrX)
    detXtX <- det(t(R) %*% R)
    reml.H0 <- -((n - p) * log(2 * pi) + (n - p) * log(sum(res^2)) + 
        log(detXtX) + (n - p) - (n - p) * log(n - p))
    #observed value of the test-statistic
    reml.obs <- 2 * logLik(m, REML = TRUE)[1]
    rlrt.obs <- max(0, reml.obs - reml.H0)
    lambda <- d$lambda
  } else {
    nonidentfixmsg <- 
      "Fixed effects structures of <mA> and <m0> not identical.
        REML-based inference not appropriate."
    if (c.m == "lme") {
      if (any(mA$fixDF$terms != m0$fixDF$terms)) 
        stop(nonidentfixmsg)
    } else {
      if (c.m == "mer") {
        if (any(mA@X != m0@X)) 
          stop(nonidentfixmsg)
      } else {
        if (c.m == "lmerMod") {
          if (any(lme4::getME(mA,"X") != lme4::getME(m0,"X")))
            stop(nonidentfixmsg)
        }
      }     
    }
    ## bug fix submitted by Andrzej Galecki 3/10/2009
    DFx <- switch(c.m, lme = anova(mA,m0)$df, 
      lmerMod = anova(mA, m0, refit = FALSE)$Df) 
    if (abs(diff(DFx)) > 1) {
      stop("Random effects not independent - covariance(s) set to 0 under H0.\n
                 exactRLRT can only test a single variance.\n")
    }
    rlrt.obs <- max(0, 2 * (logLik(mA, REML = TRUE)[1] - 
        logLik(m0, REML = TRUE)[1]))
  }
  p <- if (rlrt.obs != 0) {
    sample <- RLRTSim(X, Z, qrX=qrX, sqrt.Sigma = chol(cov2cor(Vr)), 
      lambda0 = 0, seed = seed, nsim = nsim, 
      log.grid.hi = log.grid.hi, 
      log.grid.lo = log.grid.lo, gridlength = gridlength, 
      parallel = match.arg(parallel), 
      ncpus = ncpus, cl = cl)
    if (quantile(sample, 0.9) == 0) {
      warning("Null distribution has mass ", mean(sample == 
          0), " at zero.\n")
    }
    mean(rlrt.obs < sample)
  } else {
    1
  }  
  RVAL <- list(statistic = c(RLRT = rlrt.obs), p.value = p, 
    method = paste("simulated finite sample distribution of RLRT.\n
                                (p-value based on", 
      nsim, "simulated values)"), sample=sample)
  class(RVAL) <- "htest"
  return(RVAL)
} 
