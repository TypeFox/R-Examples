# ---------------------------------------------------------------------------
# power of the TOST procedure for 2x2 crossover via simulations
#
# Author D. Labes
# ---------------------------------------------------------------------------
power.TOST.sim <- function(alpha = 0.05, logscale = TRUE, theta1, theta2, 
                           theta0, CV, n, design = "2x2", robust = FALSE, 
                           setseed = TRUE, nsims=1E5) 
{
  # Error handling
  if (missing(CV))     stop("CV must be given!")
  if (length(CV) != 1) stop("CV must be a scalar!")
  if (CV <= 0)         stop("CV must be > 0")
  if (missing(n))      stop("Number of subjects n must be given!")
  # get design characteristics
  d.no  <- .design.no(design)
  if (is.na(d.no)) stop("Design ", design, " unknown!", call. = FALSE)
  ades  <- .design.props(d.no)
  dfe   <- .design.df(ades, robust = robust)
  steps <- ades$steps
  # handle n given as scalar = ntotal or as vector = number of subjects in
  # (sequence) groups
  if (length(n) == 1) {
    if (is.finite(n)) 
      n <- nvec(n = n, grps = ades$steps)
    else
      n <- rep(Inf, times = steps)
    if (n[1] != n[length(n)]) {
      message("Unbalanced design. n(i)=", paste(n, collapse = "/"), " assumed.")
    }
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
    if (any(n < 1)) stop("All n(i) have to be >0.")
  }
  nc <- sum(1/n)
  n  <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  
  if (logscale) { # log-transformed
    if (missing(theta0))  theta0 <- 0.95
    if (any(theta0 <= 0)) stop("theta0 must be > 0.")
    if (missing(theta1))  theta1 <- 0.8
    if (missing(theta2))  theta2 <- 1/theta1
    if (theta1 <= 0 || theta1 > theta2 || length(theta1) != 1 || length(theta2) != 1)
      stop("theta1 and/or theta2 not correctly specified.")
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff <- log(theta0)
    se <- CV2se(CV) * se.fac
    mse <- CV2mse(CV)
  } else { # original scale
    if (missing(theta1))  theta1 <- -0.2
    if (missing(theta0))  theta0 <- 0.05
    if (missing(theta2))  theta2 <- -theta1
    if (theta1 > theta2 || length(theta1) != 1 || length(theta2) != 1)
      stop("theta1 and/or theta2 not correctly specified.")
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff <- theta0
    se <- CV * se.fac
    mse <- CV*CV
  }
  df <- eval(dfe)
  if (any(df < 1)) stop("n too small. Degrees of freedom <1!")

  pwr <- .power.TOST.sim(alpha = alpha, ltheta1 = ltheta1, ltheta2 = ltheta2, 
                         diffm = ldiff, mse=mse, se.fac=se.fac, df = df,
                         setseed = setseed, nsims=nsims)
  pwr
}

.power.TOST.sim <- function(alpha=0.05, ltheta1, ltheta2, diffm , mse, se.fac, 
                            df, setseed=TRUE, nsims)
{
  if(setseed) set.seed(123456)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=diffm, sd=se.fac*sqrt(mse))
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
  # the 2 one-sided test statistics
  t1 <- (pes-ltheta1)/sqrt(mses)/se.fac
  t2 <- (pes-ltheta2)/sqrt(mses)/se.fac
  # pvalues
  p1 <- 1-pt(t1, df)
  p2 <- pt(t2, df)
  #BE decision
  BE <- (p1<alpha & p2<alpha)
  
  # check against CI, debug code
#   tval <- qt(1-alpha, df)
#   hw   <- tval*sqrt(mses)*se.fac
#   lwr  <- pes - hw
#   upr  <- pes + hw
#   browser()
#   BE2  <- lwr>ltheta1 & upr<ltheta2
  # seems to give the same answer as pvalues, even if alpha>=0.5
  # will result in lwr>upr, which is not easily interpretable
  
  # number of studies found BE to number simulated
  pBE <- sum(BE)/nsims
  pBE
}