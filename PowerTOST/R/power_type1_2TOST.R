## Benjamin Lang
power.2TOST <- function(alpha = c(0.05, 0.05), logscale = TRUE, theta1, theta2,
                        theta0, CV, n, rho, design = "2x2", robust = FALSE, 
                        setseed = TRUE) {
  prob.2TOST(alpha = alpha, logscale = logscale, theta1 = theta1, 
             theta2 = theta2, theta0 = theta0, CV = CV, n = n, rho = rho,
             design = design, robust = robust, setseed = TRUE)
}

type1error.2TOST <- function(alpha = c(0.05, 0.05), logscale = TRUE, theta1, 
                             theta2, CV, n, rho, design = "2x2", robust = FALSE,
                             setseed = TRUE, details = FALSE) {
  prob.2TOST(alpha = alpha, logscale = logscale, theta1 = theta1, 
             theta2 = theta2, CV = CV, n = n, rho = rho, design = design, 
             robust = robust, setseed = TRUE, t1e = TRUE, details = details)
}

prob.2TOST <- function(alpha = c(0.05, 0.05), logscale = TRUE, theta1, theta2, 
                       theta0, CV, n, rho, design = "2x2", robust = FALSE, 
                       setseed = TRUE, t1e = FALSE, details = FALSE) {
  # Computes Power or Type I Error rate for two simultaneous TOST procedures,
  # where the two parameters of the two TOSTs are correlated with correlation
  # coefficient rho.
  #
  # Args:
  #   alpha: Vector of one-sided alpha levels for each procedure
  #   logscale: logical; if TRUE, log-transformed data is assumed
  #   theta1: Vector of lower (bio-)equivalence limits
  #   theta2: Vector of upper (bio-)equivalence limits
  #   theta0: Vector of 'true' assumed ratio of geometric means
  #   CV: Vector of coefficient of variations (use e.g. 0.3 for 30%)
  #   n: For balanced allocation to treatment/sequence groups n is the total 
  #      sample size. For unbalanced allocation n is a vector with sample size
  #      per treatment/sequence group (e.g. n = c(8, 10) for an unbalanced
  #      2x2 design with total sample size 18)
  #   rho: Correlation between the two parameters under consideration
  #   design: Character string describing the study design
  #   robust: logical; if TRUE, use robust degrees of freedom (n - #sequences)
  #   setseed: pmvt() is based on randomized quasi Monte Carlo methods.
  #            Set seed value for (pseudo) random number generator?
  #   t1e: logical; if TRUE, Type I Error will be calculated
  #   details: logical; if TRUE, return P(Type I error) for all intersection
  #            null sets, max P(Type I error) otherwise
  #
  # Returns:
  #   Type I Error if theta0 contains theta1 or theta2, Power otherwise
  #
  # Error handling
  if (length(alpha) != 2)
    stop("Two alpha values must be given!")
  if (!missing(theta0) && length(theta0) != 2)
    stop("Two theta0 values must be given!")
  if (!missing(theta1) && length(theta1) != 2)
    stop("Two theta1 values must be given!")
  if (!missing(theta2) && length(theta2) != 2)
    stop("Two theta2 values must be given!")
  if (missing(CV)) 
    stop("CV must be given!")
  if (length(CV) != 2)
    stop("Two CVs must be given!")
  if (any(CV <= 0))
    stop("CV must be > 0")
  if (missing(n)) 
    stop("Number of subjects n must be given!")
  if (missing(rho))
    stop("Correlation between the two endpoints must be given!")
  if (length(rho) != 1)
    stop("One rho must be given!")
  if (rho < -1 || rho > 1)
    stop("Correlation must be >= -1 and =< 1.")
  d.no <- .design.no(design)
  if (is.na(d.no)) 
    stop("Design ", design, " unknown!", call. = FALSE)
  ades  <- .design.props(d.no)
  dfe   <- .design.df(ades, robust = robust)
  steps <- ades$steps
  if (length(n) == 1) {
    if (is.finite(n)) 
      n <- nvec(n = n, grps = ades$steps)
    else
      n <- rep(Inf, times = steps)
    if (n[1] != n[length(n)]) {
      message("Unbalanced design. n(i)=", paste(n, collapse = "/"),
              " assumed.")
    }
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
    if (any(n < 1)) 
      stop("All n(i) have to be >0.")
  }
  nc <- sum(1/n)
  n <- sum(n)
  se.fctr <- sqrt(ades$bkni * nc)
  if (logscale) {
    if (missing(theta0))  theta0 <- c(0.95, 0.95)
    if (any(theta0 <= 0)) stop("theta0 must be > 0.")
    if (missing(theta1))  theta1 <- c(0.8, 0.8)
    if (missing(theta2))  theta2 <- 1/theta1
    if (any(theta1 <= 0) || any(theta1 > theta2))
      stop("theta1 and/or theta2 not correctly specified.")
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff <- log(theta0)
    se <- CV2se(CV) * se.fctr
  } else {
    if (missing(theta0)) theta0 <- c(0.05, 0.05)
    if (missing(theta1)) theta1 <- c(-0.2, -0.2)
    if (missing(theta2)) theta2 <- -theta1
    if (any(theta1 > theta2))
      stop("theta1 and/or theta2 not correctly specified.")
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff <- theta0
    se <- CV * se.fctr
  }
  df <- eval(dfe)
  if (any(df < 1)) 
    stop("n too small. Degrees of freedom <1!")
  if (t1e) {
    # Calculate Type I Error
    nullsets <- vector("list", 8)
    lim <- 100  # 100 instead of Inf; suffices and avoids potential
                # optimization problems when using Inf
    H_A01 <- c(-lim, log(theta1[1]))
    H_A02 <- c(log(theta2[1]), lim)
    H_C01 <- c(-lim, log(theta1[2]))
    H_C02 <- c(log(theta2[2]), lim)
    H_Aa <- c(log(theta1[1]), log(theta2[1]))
    H_Ca <- c(log(theta1[2]), log(theta2[2]))
    nullsets[[1]] <- rbind(H_A01, H_Ca)
    nullsets[[2]] <- rbind(H_A02, H_Ca)
    nullsets[[3]] <- rbind(H_Aa, H_C01)
    nullsets[[4]] <- rbind(H_Aa, H_C02)
    nullsets[[5]] <- rbind(H_A01, H_C01)
    nullsets[[6]] <- rbind(H_A01, H_C02)
    nullsets[[7]] <- rbind(H_A02, H_C01)
    nullsets[[8]] <- rbind(H_A02, H_C02)
    size.H <- function(H) {  # size of test (supremum over H)
      # Determine starting values (theta0[1], theta0[2])
      starting <- c(.getStart(H[1, ], lim), .getStart(H[2, ], lim))
      # Perform maximization over intersection nullset H via optim()
      res <- optim(par = starting, fn = .prob.2TOST, se = se, df = df, 
                   ltheta1 = ltheta1, ltheta2 = ltheta2, rho = rho, 
                   alpha = alpha, setseed = setseed, eps = 1e-05, 
                   method = "L-BFGS-B", lower = c(H[1, 1], H[2, 1]), 
                   upper = c(H[1, 2], H[2, 2]), control = list(fnscale = -1))
      if (!(res$convergence %in% c(0, 1)))
        warning("Result of maximization over nullset may not be reliable.",
                call. = FALSE) 
      # Select maximum
      res.argmax <- if (logscale) exp(res$par) else res$par
      res.max <- res$value
      c(res.max, res.argmax)
    }  # End of size.H
    # Combine results
    probs <- as.data.frame(t(vapply(nullsets, size.H, numeric(3))))
    nullsets.label <- c("H_A01 n H_Ca", "H_A02 n H_Ca", "H_Aa n H_C01",
                        "H_Aa n H_C02", "H_A01 n H_C01", "H_A01 n H_C02",
                        "H_A02 n H_C01", "H_A02 n H_C02")
    probs <- cbind(nullsets.label, probs)
    colnames(probs) <- c("Intersection null", "P(Type I Error)", "theta0 #1", 
                         "theta0 #2")
    if (details)
      return(probs)
    else
      return(max(probs[["P(Type I Error)"]]))
  } else {
    # Calculate Power
    fcn.max <- .prob.2TOST(ltheta0 = log(theta0), se = se, df = df, 
                            ltheta1 = ltheta1, ltheta2 = ltheta2, rho = rho, 
                            alpha = alpha, setseed = setseed)
    fcn.max
  }
}

.prob.2TOST <- function(ltheta0, se, df, ltheta1, ltheta2, rho, 
                        alpha = c(0.05, 0.05), setseed = TRUE, eps = 1e-04) {
  if (setseed)
    set.seed(12345)
  tval <- qt(1 - alpha, df)
  lower <- c(tval[1], -Inf, tval[2], -Inf)
  upper <- c(Inf, -tval[1], Inf, -tval[2])
  delta <- c((ltheta0[1] - ltheta1[1]) / se[1],
             (ltheta0[1] - ltheta2[1]) / se[1],
             (ltheta0[2] - ltheta1[2]) / se[2],
             (ltheta0[2] - ltheta2[2]) / se[2])
  corr <- matrix(1, ncol = 4, nrow = 4)
  corr[1:2, 3:4] <- corr[3:4, 1:2] <- rho
  prob <- pmvt(lower = lower, upper = upper, delta = delta, df = df, corr = corr,
               algorithm = GenzBretz(maxpts = 100000L, abseps = eps))
  if (attr(prob, which = "msg") != "Normal Completion")
    warning("pmvt returned message ", attr(prob, which = "msg"), call. = FALSE)
  prob[1]
}

.getStart <- function(H, lim) {
  if (H[1] <= -lim) {
    return(H[2])
  } else if (H[2] >= lim) {
    return(H[1])
  } else {
    return((H[1] + H[2]) / 2)
  }
}

# Gives the same numbers but is a bit slower
#.prob.2TOST <- function(ltheta0, se, df, ltheta1, ltheta2, rho,
#                         alpha = c(0.05, 0.05), setseed = TRUE, eps = 1e-04) {
#  if (setseed)
#    set.seed(12345)
#  tval <- qt(1 - alpha, df)
#  lower <- c(tval[1], tval[1], tval[2], tval[2])
#  upper <- c(Inf, Inf, Inf, Inf)
#  delta <- c((ltheta0[1] - ltheta1[1]) / se[1],
#             -(ltheta0[1] - ltheta2[1]) / se[1],
#             (ltheta0[2] - ltheta1[2]) / se[2],
#             -(ltheta0[2] - ltheta2[2]) / se[2])
#  corr <- diag(4)
#  corr[upper.tri(corr)] <- corr[lower.tri(corr)] <- c(-1, rho, -rho, -rho, rho, -1)
#  prob <- pmvt(lower = lower, upper = upper, delta = delta, df = df, corr = corr,
#               algorithm = GenzBretz(maxpts = 100000L, abseps = eps))
#  if (attr(prob, which = "msg") != "Normal Completion")
#    warning("pmvt returned message ", attr(prob, which = "msg"), call. = FALSE)
#  prob[1]
#}