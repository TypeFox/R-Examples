#' Likelihood Ratio Tests for simple linear mixed models
#' 
#' This function provides an exact likelihood ratio test based on simulated
#' values from the finite sample distribution for simultaneous testing of the
#' presence of the variance component and some restrictions of the fixed
#' effects in a simple linear mixed model with known correlation structure of
#' the random effect and i.i.d. errors.
#' 
#' The model under the alternative must be a linear mixed model
#' \eqn{y=X\beta+Zb+\varepsilon}{y=X*beta+Z*b+epsilon} with a \emph{single}
#' random effect \eqn{b} with known correlation structure and error terms that
#' are i.i.d. The hypothesis to be tested must be of the form \deqn{H_0:
#' \beta_{p+1-q}=\beta^0_{p+1-q},\dots,\beta_{p}=\beta^0_{p};\quad }{H0:
#' beta_1=beta0_1,..,beta_q=beta0_q, Var(b)=0}\deqn{Var(b)=0}{H0:
#' beta_1=beta0_1,..,beta_q=beta0_q, Var(b)=0} versus \deqn{H_A:\;
#' \beta_{p+1-q}\neq \beta^0_{p+1-q}\;\mbox{or}\dots }{H0: beta_1 \neq
#' beta0_1,..or..,beta_q \neq beta0_q ot
#' Var(b)>0}\deqn{\mbox{or}\;\beta_{p}\neq
#' \beta^0_{p}\;\;\mbox{or}\;Var(b)>0}{H0: beta_1 \neq beta0_1,..or..,beta_q
#' \neq beta0_q ot Var(b)>0} We use the exact finite sample distribution of the
#' likelihood ratio test statistic as derived by Crainiceanu & Ruppert (2004).
#' 
#' @param m The fitted model under the alternative; of class \code{lme},
#' \code{lmerMod} or \code{spm}
#' @param m0 The fitted model under the null hypothesis; of class \code{lm}
#' @param seed Specify a seed for \code{set.seed}
#' @param nsim Number of values to simulate
#' @param log.grid.hi Lower value of the grid on the log scale. See
#' \code{\link{exactLRT}}.
#' @param log.grid.lo Lower value of the grid on the log scale. See
#' \code{\link{exactLRT}}.
#' @param gridlength Length of the grid. See \code{\link{LRTSim}}.
#' @param parallel The type of parallel operation to be used (if any). If
#' missing, the default is "no parallelization").
#' @param ncpus integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs. Defaults to
#' 1, i.e., no parallelization.
#' @param cl An optional parallel or snow cluster for use if parallel = "snow".
#' If not supplied, a cluster on the local machine is created for the duration
#' of the call.
#' @return A list of class \code{htest} containing the following components:
#' \itemize{
#' \item \code{statistic} the observed likelihood ratio
#' \item \code{p} p-value for the observed test statistic
#' \item \code{method} a character string indicating what type of test was
#' performed and how many values were simulated to determine the critical value
#' \item \code{sample} the samples from the null distribution returned by
#' \code{\link{LRTSim}}
#' }
#' @author Fabian Scheipl, updates for \pkg{lme4.0}-compatibility by Ben Bolker
#' @seealso \code{\link{LRTSim}} for the underlying simulation algorithm;
#' \code{\link{RLRTSim}} and \code{\link{exactRLRT}} for restricted likelihood
#' based tests
#' @references Crainiceanu, C. and Ruppert, D. (2004) Likelihood ratio tests in
#' linear mixed models with one variance component, \emph{Journal of the Royal
#' Statistical Society: Series B},\bold{66},165--185.
#' @keywords htest
#' @examples
#' 
#' library(nlme);
#' data(Orthodont);
#' 
#' ##test for Sex:Age interaction and Subject-Intercept
#' mA<-lme(distance ~ Sex * I(age - 11), random = ~ 1| Subject,
#'   data = Orthodont, method = "ML")
#' m0<-lm(distance ~ Sex + I(age - 11), data = Orthodont)
#' summary(mA)
#' summary(m0)
#' exactLRT(m = mA, m0 = m0)
#' 
#' @export exactLRT
#' @importFrom stats coefficients
`exactLRT` <-
  function(m, m0, seed = NA, nsim = 10000, 
    log.grid.hi = 8, log.grid.lo = -10, gridlength = 200,
    parallel = c("no", "multicore", "snow"), 
    ncpus = 1L, cl = NULL) 
  {
    if (class(m0) != "lm") 
      stop("m0 not an lm-object. \n")
    if (class(m) == "spm") {
      m <- m$fit
      class(m) <- "lme"
    }
    if (class(m) %in% c("amer", "mer"))
      stop("Models fit with package <amer> or versions of <lme4> below 1.0 are no longer supported.")
    if (!((c.m <- class(m)) %in% c("lme", "lmerMod", "merModLmerTest"))) 
      stop("Invalid <m> specified. \n")
    
    
    d <- switch(c.m, lme = extract.lmeDesign(m), 
      lmerMod=extract.lmerModDesign(m))
    if(length(d$lambda) != 1 || d$k != 1) 
      stop("multiple random effects in model - 
                 exactLRT needs <m> with only a single random effect.")
    X <- d$X
    Z <- d$Z
    y <- d$y
    Vr <- d$Vr
    K <- NCOL(Z)
    n <- NROW(X)
    p <- NCOL(X)
    q <- p - length(coefficients(m0)[!is.na(coefficients(m0))])
    if (n != length(m0$fitted)) 
      stop("different data under the null and alternative. \n")
    if (q < 0) 
      stop("m0 not nested in m. \n")
    if (n - p - K < 1) 
      stop("No. of effects greater than no. of observations. Reduce model complexity.\n")
    if (q == 0) 
      message("No restrictions on fixed effects. REML-based inference preferable.")
    method <- switch(c.m, lme = m$method, 
      lmerMod=ifelse(lme4::isREML(m), "REML", "ML"))
    if (method != "ML") {
      message("Using likelihood evaluated at REML estimators.")
      message("Please refit model with method=\"ML\" for exact results.")
    }
    #observed value of the LRT
    lrt.obs <- max(0, 2 * logLik(m, REML = FALSE)[1] - 2 * logLik(m0, 
      REML = FALSE)[1])
    sample <- LRTSim(X, Z, q, sqrt.Sigma = chol(cov2cor(Vr)), 
      seed = seed, nsim = nsim, log.grid.hi = log.grid.hi, 
      log.grid.lo = log.grid.lo, gridlength = gridlength,
      parallel = match.arg(parallel), 
      ncpus = ncpus, cl = cl)
    if (quantile(sample, 0.9) == 0) {
      warning("Null distribution has mass ", mean(sample == 
          0), " at zero.\n")
    }
    p <- mean(lrt.obs < sample)
    RVAL <- list(statistic = c(LRT = lrt.obs), 
      p.value = p, 
      method = paste("simulated finite sample distribution of LRT. (p-value based on", 
        nsim, "simulated values)"), sample=sample)
    class(RVAL) <- "htest"
    return(RVAL)
  } 

