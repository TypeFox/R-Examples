#-----------------------------------------------------------------------
# Sample size for partial and full replicate design and scaled ABE
# via simulated (empirical) power. Alpha is adjusted to maintain the
# empiric TIE <= nominal alpha.
#
# Author: Helmut Schuetz
#-----------------------------------------------------------------------
sampleN.scABEL.ad <- function(alpha = 0.05, targetpower = 0.8, theta0,
                              theta1, theta2, CV = 0.3,
                              design = c("2x3x3", "2x2x4", "2x2x3"),
                              regulator = c("EMA", "ANVISA"),
                              nstart = NA, nsims = 1e6, imax=100,
                              tol, print = TRUE, details = FALSE,
                              alpha.pre = 0.05, setseed = TRUE)
{
  ## Arguments:
  ##   alpha       Nominal alpha (in BE generally fixed to 0.05).
  ##               Lower value only if needed (e.g. to correct for
  ##               multiplicity).
  ##   targetpower Desired power.
  ##   theta0      Expected GMR. Defaults to 0.9 - different from the
  ##               default 0.95 in sampleN.scABEL()!
  ##   theta1      Lower margin. Defaults to 0.8.
  ##   theta2      Upper margin. Defaults to 1/theta1.
  ##   CV          Intra-subject CV(s) obtained in a replicate design.
  ##               (ratio, /not/ percent).
  ##               If given as a scalar, the CV of R.
  ##               If given as a vector, CV[1] /must/ be the CV of T and
  ##               CV[2] the CV of R. Important!
  ##   design      "2x2x4", "2x2x3", "2x3x3"
  ##   regulator   "EMA" or "ANVISA"
  ##               Cave: ANVISA's requirements are unofficial and might
  ##               require extreme adjustment close to CVwR 40%.
  ##   nstart      If given, the starting sample size.
  ##   nsims       Simulations for the TIE. Should not be <1e6.
  ##   imax        max. number of steps in sample size search
  ##   tol         desired accuracy (convergence tolerance)
  ##               defaults to 1e-6 for EMA and 1e-7 for ANVISA
  ##   print       Boolean. If FALSE, returns a data.frame of results.
  ##   details     Boolean (intermediates, runtime, number of sim's).
  ##   alpha.pre   Pre-specified level.
  ##   setseed     Boolean (default TRUE uses set.seed(123456)).
  ## Returns:
  ##   n           Sample size which maintains the TIE for the
  ##               adjusted (or pre-specified) alpha.
  ##   alpha.adj   Iteratively adjusted alpha.
  ##   pwr         Achieved power.
  ##   TIE         Empiric Type I Error (aka rejection rate).
  ## Algo:
  ##   1. Estimate the TIE for the /unadjusted/ alpha (alpha) or
  ##      the /pre-specified/ alpha (alpha.pre).
  ##   2. If no inflation of TIE, stop. Othewise, continue with 3 - 5.
  ##   3. Iteratively adjust alpha to preserve the consumer risk.
  ##   4. Get a new sample size for /this/ alpha (might be higher) and
  ##      estimate the TIE.
  ##   5. Increase the sample size and repeat steps 3 & 4 until the
  ##      target power is reached.
  ######################################################################
  ## Tested on Win 7 Pro SP1 64bit
  ##   R 3.2.4 Revised 64bit (2016-03-16), PowerTOST 1.3-4 (2016-03-10)
  ######################################################################
  env <- as.character(Sys.info()[1]) # get info about the OS
  if ((env == "Windows") || (env == "Darwin")) flushable <- TRUE
    else flushable <- FALSE # supress flushing on other OS's
  # acceptance range defaults
  if (missing(theta1) && missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 = 1/theta1
  # check theta0
  if (missing(theta0)) theta0 <- 0.9
  if (theta0 < theta1 || theta0 > theta2)
    stop("theta0 must be within [theta1, theta2]")
  # check regulator arg
  regulator <- toupper(regulator)
  regulator <- match.arg(regulator)
  if (length(nstart) == 2) nstart <- sum(nstart)
  design <- match.arg(design)
  CVwT <- CV[1]
  if (length(CV) == 2) CVwR <- CV[2] else CVwR <- CVwT
  if (!is.na(nstart) &&
    ((regulator == "EMA" && nstart < 12) ||
     (regulator == "ANVISA" && nstart < 24)))
      warning("Requested sample size below regulatory minimum.")
  if (!is.na(targetpower) && (targetpower < 0 || targetpower >= 1))
    stop("targetpower must be within 0 <= 1.")
  if (alpha.pre > alpha) {
    warning(paste0("alpha.pre > alpha doesn't make sense.",
                   "\nalpha.pre was set to alpha."))
    alpha.pre <- alpha
  }
  seqs <- as.numeric(substr(design, 3, 3))  # subjects / sequence
  sig  <- binom.test(x = round(alpha*nsims, 0), n = nsims,
                     alternative = "less",
                     conf.level = 1 - alpha)$conf.int[2]
  method <- "ABE"
  if ((regulator == "EMA" && CVwR > 0.3) ||
      (regulator == "ANVISA" && CVwR > 0.4)) method <- "ABEL"
  # define the data.frame of rseults
  res <- data.frame(design = design, regulator = regulator,
                    method = method, theta0 = theta0, CVwT = CVwT,
                    CVwR = CVwR, alpha = alpha, alpha.pre = alpha.pre,
                    alpha.adj = NA, TIE = NA, n = NA,
                    targetpower = targetpower, power = NA)
  names(res) <- c("Design", "Regulator", "Method", "theta0", "CVwT",
                  "CVwR", "alpha", "alpha.pre", "adj. alpha", "TIE",
                  "Sample size", "Target power", "Achieved power")
  limits <- as.numeric(scABEL(CV = CVwR, regulator = regulator))
  U <- limits[2] # Simulate at the upper (expanded) limit. For CVwR
                 # 30% that's 1.25. Due to the symmetry simulations
                 # at the lower limit (0.8) would work as well.
  if (is.na(alpha.pre) || (alpha.pre != alpha)) {
    al <- alpha.pre # If pre-specified, use alpha.pre.
  } else {
    al <- alpha     # If not, use alpha (commonly 0.05).
  }
  designs <- c("2x2x4", "2x2x3", "2x3x3")
  type    <- c("RTRT|TRTR", "RTR|TRT", "RRT|RTR|TRR") # clear words
  if (print) { # Show input to keep the spirits of the user high.
    cat("\n+++++++++++ scaled (widened) ABEL +++++++++++\n")
    cat("            Sample size estimation\n")
    cat("        for iteratively adjusted alpha\n")
    cat("---------------------------------------------\n")
    cat("Study design: ")
    cat(paste0(design, " (", type[match(design, designs)], ")\n"))
    cat("log-transformed data (multiplicative model)\n")
    cat(formatC(nsims, format = "d", big.mark = ",", decimal.mark = "."),
        "studies in each iteration simulated.\n\n")
    txt <- paste0("Expected CVwR ", sprintf("%.4g", CVwR))
    if (length(CV) == 2) {
      txt <- paste0(txt, ", CVwT ", sprintf("%.4g", CVwT), "\n")
    } else {
      txt <- paste0(txt, "\n")
    }
    cat(txt)
    txt <- paste0("Nominal alpha      : ", signif(alpha, 5))
    if (!is.na(alpha.pre) && (alpha.pre != alpha)) {
      txt <- paste0(txt, ", pre-specified alpha ", alpha.pre, "\n")
    } else {
      txt <- paste(txt, "\n")
    }
    cat(txt)
    cat("Null (true) ratio  :", sprintf("%.3f", theta0), "\n")
    cat("Target power       :", sprintf("%.3g", targetpower), "\n")
    cat(paste0("Regulatory settings: ", regulator, " (", method, ")\n"))
    if (regulator == "EMA") {
      if (CVwR <= 0.3) {
        cat("Switching CVwR     : 0.30",
            "\nBE limits          : 0.8000...1.2500\n")
      } else {
        cat(paste("Switching CVwR     : 30%",
                  "\nRegulatory constant: 0.760\n"))
        cat(sprintf("%s    : %.4f%s%.4f%s", "Expanded limits",
                    limits[1], "...", limits[2], "\n"))
        cat("Upper scaling cap  : CVwR 0.5\n")
        cat("PE constraints     : 0.8000...1.2500\n")
      }
    } else {
      if (CVwR <= 0.4) {
        cat("Switching CVwR     : 0.40 (unofficial)",
            "\nBE limits          : 0.8000...1.2500\n")
      } else {
        cat("Switching CVwR     : 0.40 (unofficial)\nRegulatory constant:",
            signif(log(1.25)/CV2se(0.3), 7),
            "(assumed; no official guidance)\n")
        cat(sprintf("%s    : %.4f%s%.4f%s", "Expanded limits",
                    limits[1], "...", limits[2], "\n"))
        cat("Upper scaling cap  : CVwR 0.50\n")
        cat("PE constraints     : 0.8000...1.2500\n")
      }
    }
    if (flushable) flush.console()
  }
  if (details) ptm <- proc.time()
  no <- 0 # Simulation counter.
  if (is.na(nstart)) { # If sample size is not given, estimate one.
    unadj.n  <- sampleN.scABEL(alpha = al, targetpower = targetpower,
                               theta0 = theta0, CV = CV, design = design,
                               regulator = regulator, imax=imax,
                               print = FALSE, details = FALSE, nsims = 1e5,
                               setseed = setseed)[["Sample size"]]
    no <- 1e5
  } else {             # Start with the specified sample size.
    unadj.n <- nstart
  }
  # Get results for the sample size.
  x <- scABEL.ad(alpha = alpha, theta0 = theta0, CV = CV,
                 design = design, regulator = regulator, n = unadj.n,
                 nsims = nsims, imax=imax, print = FALSE, details = FALSE,
                 alpha.pre = alpha.pre, setseed = setseed)
  alpha.adj <- x[["alpha.adj"]]
  if (is.na(alpha.adj)) { # No adjustment was necessary:
    if (alpha.pre != alpha) {
      alpha.adj <- alpha.pre # If pre-specified, use alpha.pre.
    } else {
      alpha.adj <- alpha     # If not, use alpha (commonly 0.05).
    }
  }
  TIE.unadj <- x[["TIE.unadj"]]
  pwr.unadj <- x[["pwr.unadj"]]
  TIE.adj   <- x[["TIE.adj"]]
  pwr.adj   <- x[["pwr.adj"]]
  no        <- no + x[["sims"]]
  if (print) {
    if (alpha.pre == alpha)
      al.txt <- "nomin. alpha:" else al.txt <- " spec. alpha:"
    if (details) {
      cat(sprintf("%s %3d, %s %.4f %s %.4f%s %.4f%s", "\nn", unadj.n,
                  al.txt, al, "(power", pwr.unadj, "), TIE:", TIE.unadj,
                  "\n"))
    }
  }
  step.1 <- FALSE # Check conditions for stopping below:
  if (TIE.unadj <= alpha && pwr.unadj > targetpower) step.1 <- TRUE
  if (!is.na(TIE.adj)) {
    if (TIE.adj <= alpha && pwr.adj > targetpower) step.1 <- TRUE
  }
  # browser()
  if (step.1 && is.na(TIE.adj)) { # Stop: Nothing to do.
    if (!details && print) { # only if we don't have this info already
      cat(sprintf("%s %3d, %s %.4f %s %.4f%s %.4f%s", "\nn", unadj.n,
                  al.txt, al, "(power", pwr.unadj, "), TIE:",
                  TIE.unadj, "\n"))
    }
    if (print) {
      cat("No inflation of the TIE expected; ")
      if (alpha.pre != alpha) {
        cat("the chosen pre-specified alpha is justified.\n")
      } else {
        cat("hence, no adjustment of alpha required.\n")
      }
    }
    res[["TIE"]]            <- signif(TIE.unadj, 5)
    res[["Sample size"]]    <- unadj.n
    res[["Achieved power"]] <- signif(pwr.unadj, 5)
    return(invisible(res))
  }
  if (print && details && (alpha.adj != alpha)) { # Some information.
    cat("\nSample size search and iteratively adjusting alpha")
    cat(sprintf("%s %3d, %s ", "\nn", unadj.n, "  adj. alpha:"))
    if (alpha.adj >= 0.01) {    # General case (alphas <0.025 are rare).
      cat(sprintf("%.5f %s %.4f%s %.2f%%%s", alpha.adj, "(power", pwr.adj,
                  "), rel. impact on power:",
                  100*(pwr.adj - pwr.unadj)/pwr.unadj, "\n"))
    } else {                # Sometimes necessary for ANVISA...
      cat(paste0(signif(alpha.adj, 5),
          sprintf(" %s %.4f%s %.2f%%%s", "(power", pwr.adj,
                  "), relative impact on power:",
                  100*(pwr.adj - pwr.unadj)/pwr.unadj, "\n")))
    }
    if (flushable) flush.console()
  }
  # Increase the sample size /and/ adjust alpha until achieved
  # power is at least the target power and the TIE does not
  # exceed (nominal) alpha.
  pwr <- iter <- 0
  while (pwr < targetpower) {
    if (iter == 0) { # Get the sample size for the (first!) adjusted
                     # alpha obtained from scABEL.ad(...) above.
                     # Faster than scABEL.ad() because only 1e5 sim's.
      n.new <- sampleN.scABEL(alpha = alpha.adj, CV = CV, theta0 = theta0,
                              targetpower = targetpower, design = design,
                              regulator = regulator, imax=imax, print = FALSE,
                              details = FALSE,
                              setseed = setseed)[["Sample size"]]
      no <- no + 1e5
      step.1 <- TRUE
    } else {         # In later iterations use scABEL.ad().
      if (step.1) {             # Prevents overshooting power in
        n.new <- n.new - 2*seqs # the 1st step, e.g, aim /lower/
        step.1 <- FALSE         # to be on the safe side!
      } else {
        n.new <- n.new + seqs   # Increase n in further steps.
      }
    }
    if (alpha.adj != alpha.pre) { # Adjust alpha (general case).
      x  <- scABEL.ad(alpha = alpha, regulator = regulator, design = design,
                      CV = CV, n = n.new, theta0 = theta0, imax=imax,
                      tol = tol, print = FALSE, details = FALSE,
                      nsims = nsims, setseed = setseed)
    } else {                      # Do /not/ adjust pre-specified alpha!
      x  <- scABEL.ad(regulator = regulator, design = design, CV = CV,
                      n = n.new, theta0 = theta0, imax=imax, tol = tol,
                      print = FALSE, details = FALSE, nsims = nsims,
                      setseed = setseed, alpha.pre = alpha.adj)
    }
    no <- no + x$sims
    if (is.na(x[["alpha.adj"]])) { # No adjustment was necessary!
      pwr       <- x[["pwr.unadj"]]
    } else {
      pwr       <- x[["pwr.adj"]]
      alpha.adj <- x[["alpha.adj"]]
      TIE       <- x[["TIE.adj"]]
    }
    iter <- iter + 1
    # browser()
    if (pwr < targetpower && iter >= 1) { # Show intermediate steps.
      if (print && details) {
        if (alpha.adj >= 0.01) { # Nice format (EMA)
          cat(sprintf("%s %3d, %s %.5f %s %.4f%s", "n", n.new,
                      "  adj. alpha:", alpha.adj, "(power", pwr, ")\n"))
        } else {                 # Sometimes needed for ANVISA.
          cat(sprintf("%s %3d, %s ", "n", n.new, "  adj. alpha:"))
          cat(paste0(signif(alpha.adj, 5),
              sprintf(" %s %.5f%s", "(power", pwr, ")\n")))
        }
        if (flushable) flush.console() # advance console output.
      }
    }
  }
  if (details) run.time <- proc.time() - ptm
  if (print) {
    cat(sprintf("%s %3d, %s ", "n", n.new, "  adj. alpha:"))
    if (alpha.adj >= 0.01) { # As above. EMA
      cat(sprintf("%.5f %s %.4f%s %.5f%s", alpha.adj, "(power", pwr,
                  "), TIE:", TIE, "\n"))
    } else {                 # As above. ANVISA
      cat(signif(alpha.adj, 5),
          sprintf("%s %.4f%s %.5f%s", "(power", pwr,
                  "), TIE:", TIE, "\n"))
    }
    if (details) {
      cat("Compared to nominal alpha's sample size increase of",
          sprintf("%.1f%%", 100*(n.new - unadj.n)/unadj.n),
          "(~study costs).\n\n")
    } else {
      cat("\n\n")
    }
  }
  if (print && details) {
    cat("Runtime    :", signif(run.time[3], 3), "seconds",
        "\nSimulations:", formatC(no, format = "d", big.mark = ",",
                                  decimal.mark = "."), "\n\n")
  }
  if (TIE > sig) { # Happens rarely for ANVISA (only).
    warning(paste0("Algorithm failed. ",
                   "Try to restart with at least 'nstart = ",
                   n.new + seqs, "'."))
  }
  res[["adj. alpha"]]     <- signif(alpha.adj, 5)
  res[["TIE"]]            <- signif(TIE, 5)
  res[["Sample size"]]    <- n.new
  res[["Achieved power"]] <- signif(pwr, 5)
  if (print || details) {
    return(invisible(res))
  } else {
    return(res)
  }
}
# Examples
#   sampleN.scABEL.ad(regulator="EMA", design="2x2x4", CV=0.3, theta0=0.9, targetpower=0.8, details=TRUE)
# should return:
#   +++++++++++ scaled (widened) ABEL +++++++++++
#              Sample size estimation
#           for iteratively adjusted alpha
#   ---------------------------------------------
#   Study design: 2x2x4 (RTRT|TRTR)
#   log-transformed data (multiplicative model)
#   1,000,000 studies in each iteration simulated.
#
#   Expected CVwR 0.3
#   Nominal alpha      : 0.05
#   Significance limit : 0.05036
#   Null (true) ratio  : 0.900
#   Regulatory settings: EMA (ABE)
#   Switching CVwR     : 0.30
#   BE limits          : 0.8000...1.2500
#
#   n  34, nomin. alpha: 0.0500 (power 0.8028), TIE: 0.0816
#
#   Sample size search and iteratively adjusting alpha
#   n  34,   adj. alpha: 0.0286 (power 0.7251), rel. impact on power: -9.68%
#   n  42,   adj. alpha: 0.0283 (power 0.8022), TIE: 0.0500
#   Compared to nominal alpha's sample size increase of 23.5% (~study costs).
#
#   Runtime    : 11.6 seconds
#   Simulations: 10,400,000
#
#   x <- sampleN.scABEL.ad(regulator="EMA", design="2x2x3", CV=c(0.35, 0.40), nstart=42, theta0=0.9, targetpower=0.8, details=FALSE, print=FALSE)
#   Show the results:
#   print(x, row.names=FALSE)
#    Design Regulator Method theta0 CVwT CVwR alpha alpha.pre adj. alpha      TIE Sample size Target power Achieved power
#     2x2x3       EMA   ABEL    0.9 0.35  0.4  0.05      0.05   0.038992 0.050001          46          0.8        0.81531
#   Show the sample size only: x[["Sample size"]]
#   [1] 46
