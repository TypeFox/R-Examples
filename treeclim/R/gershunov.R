##' Gershunov test for spurious low-frequency modulations
##' 
##' This function provides a test to decide whether low-frequency
##' modulations in the relationship between climate and tree-growth
##' are significantly stronger or weaker than could be expected by
##' chance.
##' @details This function is a multivariate extension of the test for
##' spurious low-frequency modulations for moving correlations of time
##' series as proposed by Gershunov et al. (2001). In short, 1000
##' simulations of random data sets are generated, where the climate
##' data is simulated as Gaussian noise, and the tree-data as linear
##' combinations of the climate parameters using the original
##' coefficients of the correlation function, and an error component
##' with a variance equal to the variance unexplained by the
##' individual parameters.
##'   
##' For each iteration, a moving correlation function is calculated
##' with exactly the same settings as the original model. The standard
##' deviation over the individual windows for each parameter is then
##' compared to the bootstrapped distribution of the standard
##' deviation of the simulated data to test for significantly higher
##' or lower low-frequency modulations.
##' @param x an object of class '"tc_dcc"' as returned from a call to
##'     \code{\link{dcc}} with moving correlations enabled
##' @param boot \code{logical} shall the individual correlation be
##'     bootstrapped?  (see details)
##' @param sb \code{logical} shall a status bar be drawn?
##' @return a \code{data.frame} with p values for the testing the null
##'     hypothesis that the low-frequency modulation of the
##'     correlations of the variables with tree-growth can be
##'     considered as noise.
##' @references Gershunov, A., N. Schneider, and
##'     T. Barnett. 2001. Low-frequency modulation of the ENSO-Indian
##'     Monsoon rainfall relationship: Signal or noise? Journal of
##'     Climate 14:2486-2492.
##' @examples
##' \dontrun{
##' dc_cor <- dcc(muc_spruce, muc_clim, 3:9, method = "cor", moving = TRUE)
##' g_test(dc_cor)
##' }
##' @keywords test
##' @export
g_test <- function(x, boot = FALSE, sb = TRUE) {
  if (!any(class(x) != "tc_dcc"))
    stop("Please provide output of function `dcc`.")

  if (is.null(x$call$moving) || x$call$moving == FALSE)
    stop("Gershunov test can only be computed for moving correlations.")
  
  if (length(pmatch(x$call$method, "correlation")) == 0)
      stop("Gershunov test is currently only implemented for running correlation functions, not for response function.")

  ## number of bootstrap resamples
  niter <- 500

  ## get parameters for moving correlation function from call
  .win_size <- ifelse(is.null(x$call$win_size), 25, x$call$win_size)
  .win_offset <- ifelse(is.null(x$call$win_offset), 1, x$call$win_offset)
  .start_last <- ifelse(is.null(x$call$start_last), TRUE, x$call$start_last)
  .boot <- ifelse(is.null(x$call$boot), "stationary", x$call$boot)
  
  ## calculate null model with full data set to get coefficients
  c0 <- tc_correlation(x$truncated$tree, x$design, ci = 0.05,
                       boot = .boot)$result$coef
  n <- length(c0)
  m <- length(x$truncated$tree)
  
  ## get unexplained variance by model
  prediction <- rowSums(t(t(scale(x$design$aggregate)) * c0))
  unex <- 1 - var(prediction)
  if (unex < 0) {
    unex <- 0
  }
  
  ## overwrite .boot, if set to FALSE for g-test (default)
  if (!boot)
    .boot <- "none"
  
  ## if boot is TRUE, warn because of long calculation, and give the
  ## option to cancel
  if (boot) {
    cat("Checking duration...\n")
    dur <- system.time({
      tc_mfunc(x$truncated$tree, x$design,
               method = "correlation",
               start_last = .start_last,
               win_size = .win_size,
               win_offset = .win_offset,
               boot = .boot,
               sb = FALSE,
               ci = 0.05)
    })
    dur1000 <- ceiling(dur[3] * niter / 60)
    cat("Running this test with bootstrapping on the individual correlations enabled will take around", dur1000, "minutes.\n")
    ans <- readline("Do you really want to run this? [Y/n]\n")
    if (ans != "Y")
      stop("Execution interrupted by user.")
    cat("Ok then. Maybe time to grab a sandwich?\n")
  } else {
    cat("Performing test without bootstrapping individual correlations.\n")
  }
  
  ## test tree-ring data and the columns of the design matrix for normality
  pn_tree <- shapiro.test(x$truncated$tree)$p.value > 0.05
  if (!pn_tree)
    warning("Tree-ring data is not normally distributed. Gershunov-Test might yield meaningless results.")
  
  pn_design <- apply(x$design$aggregate, 2,
                     function(x) shapiro.test(x)$p.value) > 0.05
  if (!all(pn_design))
    warning("Not all climate variables are normally distributed. Gershunov-Test might yield meaningless results.")
  
  ## get sd of parameters for original run
  sd0 <- apply(x$coef$coef, 1, sd)
  
  ## for each c0, simulate niter pairs of random time series
  sds <- matrix(NA, ncol = niter, nrow = n)
  
  if (sb)                            # initialize status bar (if TRUE)
      mpb <- txtProgressBar(min = 1,  max = niter, style = 3)

  for (i in 1:niter) {
    ## model tree-ring series (gamma2) as linear combination of
    ## climate (gamma1) + error term representing the variance
      ## unexplained by original model
    gamma1 <- matrix(rnorm(n * m, 0, 1), nrow = m)
    rownames(gamma1) <- rownames(x$design$aggregate)
    gamma2 <- rowSums(t(t(gamma1) * c0)) + rnorm(m, 0, sd = sqrt(unex))
      
    gamma1l <- list(aggregate = gamma1,
                   names = x$design$names,
                   pretty_names = x$design$pretty_names)
    
    ## throw data into the _same_ moving correlation function
    c1 <- tc_mfunc(gamma2, gamma1l,
                   method = "correlation",
                   start_last = .start_last,
                   win_size = .win_size,
                   win_offset = .win_offset,
                   boot = .boot,
                   sb = FALSE,
                   ci = 0.05)$result$coef
    
    sds[,i] <- apply(c1, 1, sd)
    
    if (sb)                             # update status bar (if TRUE)
      setTxtProgressBar(mpb, i)
  }  
  
  ## calculate p values from cumulative density function
  ps <- numeric(n)
  for (i in 1:n) {
      ps[i] <- ecdf(sds[i,])(sd0[i])
  }
  
  ## prepare output
  out <- data.frame(
    varname = rev(x$design$pretty_names$varname),
    month = rev(x$design$pretty_names$month_label),
    p = rev(1 - ps)
    )
  
  if (sb)                               # close status bar (if TRUE)
    close(mpb)
  
  out
}
