print.summary.limitmeta <- function(x,
                                    backtransf = x$backtransf,
                                    digits = max(3, .Options$digits - 3),
                                    header = TRUE, ...){
  
  meta:::chkclass(x, "summary.limitmeta")
  
  
  sm <- x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summary.limitmeta"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  
  
  sm.lab <- sm
  ##
  if (backtransf){
    if (sm == "ZCOR")
      sm.lab <- "COR"
    if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT", "PLN"))
      sm.lab <- "proportion"
  }
  else 
    if (meta:::is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
  
  
  TEs <- c(x$TE.adjust, x$TE.random)
  lower <- c(x$lower.adjust, x$lower.random)
  upper <- c(x$upper.adjust, x$upper.random)


  if (backtransf){
    ##
    npft.ma <- 1 / mean(1 / x$x$n)
    ##
    TEs   <- meta:::backtransf(TEs, sm, "mean",
                               npft.ma, warn = TRUE)
    lower <- meta:::backtransf(lower, sm, "lower",
                               npft.ma, warn = TRUE)
    upper <- meta:::backtransf(upper, sm, "upper",
                               npft.ma, warn = TRUE)
  }
  ##
  TEs   <- round(TEs, digits)
  lower <- round(lower, digits)
  upper <- round(upper, digits)
  ##
  TEs   <- format(TEs, trim = TRUE)
  lower <- format(lower, trim = TRUE)
  upper <- format(upper, trim = TRUE)
  ##
  pvals <- c(x$pval.adjust, x$pval.random)
  zvals <- format(round(c(x$zval.adjust, x$zval.random), digits))
  
  
  imeth <- charmatch(tolower(x$method.adjust),
                     c("beta0", "betalim", "mulim"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("Argument 'method.adjust' should be \"beta0\", \"betalim\", or \"mulim\"")
  ##
  method.adjust.detail <- c("- expectation (beta0)",
                            "- including bias parameter (beta-lim)",
                            "- excluding bias parameter (mu-lim)")[imeth]
  
  
  ci.lab <- paste(round(100 * x$level.comb, 1),
                  "%-CI", sep="")
  
  
  res <- cbind(c("Adjusted estimate",
                 "Unadjusted estimate"),
               ifelse(TEs == "NA", "", TEs),
               meta:::p.ci(lower, upper),
               zvals,
               meta:::format.p(pvals)
               )
  ##
  dimnames(res) <- list(rep("", dim(res)[[1]]),
                        c("Random effects model",
                          sm.lab, ci.lab, "z", "pval"))
  
  
  if (header)
    meta:::crtitle(x)
  
  
  cat("Result of limit meta-analysis:\n\n")
  ##
  prmatrix(res, quote = FALSE, right = TRUE)
  
  
  if (!is.na(x$tau)){
    ##
    I2res <- meta:::isquared(x$Q, x$k - 1, x$x$level.comb)
    I2 <- I2res$TE
    lowI2 <- I2res$lower
    uppI2 <- I2res$upper
    ##
    cat(paste("\nQuantifying heterogeneity:\n",
              if (x$tau^2 > 0 & x$tau^2 < 0.0001)
              "tau^2 < 0.0001"
              else
              paste("tau^2 = ",
                    ifelse(x$tau == 0,
                           "0",
                           format(round(x$tau^2, 4), 4, nsmall = 4, scientific = FALSE)),
                    sep="")
              ,
              paste("; I^2 = ", round(100 * I2, 1), "%",
                    ifelse(x$k > 2,
                           paste(" ",
                                 meta:::p.ci(paste(round(100 * lowI2, 1), "%", sep = ""),
                                             paste(round(100 * uppI2, 1), "%", sep = "")),
                                 sep = ""),
                           ""),
                    "; G^2 = ", round(100 * x$G.squared, 1), "%",
                    sep = ""),
              "\n", sep = ""))
  }
  
  
  Qd <- function(Q, df)
    data.frame(Q = round(Q, 2), df = df,
               p = 1 - pchisq(Q, df = df))
  ##
  Qds <- rbind(Qd(x$Q, x$k - 1),
               Qd(x$Q.small, 1),
               Qd(x$Q.resid, x$k - 2))
  Qds$Q  <- format(Qds$Q)
  Qds$df <- format(Qds$df)
  Qds$p  <- meta:::format.p(Qds$p)
  ##
  Qd1 <- Qds[1, ]
  Qd2 <- Qds[2, ]
  Qd3 <- Qds[3, ]
  ##
  dimnames(Qd1) <- list("", c("Q", "d.f.", "p-value"))
  cat("\nTest of heterogeneity:\n")
  prmatrix(Qd1, quote = FALSE, right = TRUE, ...)
  ##
  dimnames(Qd2) <- list("", c("Q-Q'", "d.f.", "p-value"))
  cat("\nTest of small-study effects:\n")
  prmatrix(Qd2, quote = FALSE, right = TRUE, ...)
  ##
  dimnames(Qd3) <- list("", c("Q'", "d.f.", "p-value"))
  cat("\nTest of residual heterogeneity beyond small-study effects:\n")
  prmatrix(Qd3, quote = FALSE, right = TRUE, ...)
  ##
  cat("\nDetails on adjustment method:\n",
      method.adjust.detail,
      "\n", sep = "")
  
  
  invisible(NULL)
}
