##===================================================================
## Retrieves colnames for confidence or prediction bands
## from a dataframe given in 'x'. The levels are given in
## the colnames
##
## NEW IN VERSION  1.3-5
##===================================================================
predNames <- function(x, prefix = c("L", "U")) {
  cn <- names(x)
  mn <- character(0)
  pct <- character(0)
  type <- character(0)
  ## find the narrowest confint
  if ("L" %in% prefix) {
    ##pct.L <- grep("L.[0-9]*$", cn, value = TRUE)
    pct.L <-  grep("^L.[0-9]*\\.?[0-9]*$", cn, value = TRUE)
    if (ln <- length(pct.L)) {
      mn <- c(mn, pct.L)
      type <- c(type, rep("L", ln))
      pct <- c(pct, gsub("L.", "", pct.L))
    }
  }
  if ("U" %in% prefix) {
    pct.U <- grep("^U.[[:digit:]]*\\.?[[:digit:]]*$", cn, value = TRUE)
    if (ln <- length(pct.U)) {
      mn <- c(mn, pct.U)
      type <- c(type, rep("U", ln))
      pct <- c(pct, gsub("U.", "", pct.U))
    }
  }
  list(names = mn,
       type = type,
       pct = as.numeric(pct))
}

##===================================================================
## Compute limits for levels using confidence or prediction bands
## as well as historical data.
##
## NEW IN VERSION 1.4-0
##
##===================================================================
rangeLev.Renouv <- function(x,
                            show.MAX = TRUE,
                            show.OTS = TRUE,
                            Tlim = NULL) {
  pn <- predNames(x$ret.lev, prefix = c("L", "U"))
  pnn <- c("quant", pn$names)
  ## cat("pn = ", pnn, "\n")
  Tlim <- range(x$pred$period, Tlim)
  ## cat("x$pred = \n"); print(x$pred)
  ## cat("Tlim = \n"); print(Tlim)
  if (!is.null(Tlim)) {
    period <- x$ret.lev$period
    ind <- (period >= Tlim[1]) & (period <= Tlim[2])
    r <- range(as.matrix(x$ret.lev[ind, pnn]), x$x.OT, na.rm = TRUE)
  } else {
    r <- range(as.matrix(x$ret.lev[ , pnn]), x$x.OT, na.rm = TRUE)
  }
  if (x$history.MAX$flag && show.MAX) {
    x1 <- unlist(x$history.MAX$data)
    if (length(x1)) r <- range(r, range(x1), na.rm = TRUE)
  } 
  if (x$history.OTS$flag && show.OTS) {
    x1 <- unlist(x$history.OTS$threshold)
    if (length(x1)) r <- range(r, range(x1), na.rm = TRUE)
    x1 <- unlist(x$history.OTS$data)
    if (length(x1)) r <- range(r, range(x1), na.rm = TRUE)
  }
  r 
}

##=====================================================================
## round quantiles and confidence limits using a 
## suitable number of digits
##=====================================================================
roundPred <- function(pred, dig.quant = NA) {
  cn <- colnames(pred)
  ## find the narrowest confint
  pct.L <- grep("L.[0-9]", cn, value = TRUE)
  pct.U <- grep("U.[0-9]", cn, value = TRUE)
  if (is.na(dig.quant) || length(dig.quant) == 0) {
    disp <- FALSE
    if (length(pct.L)) {
      i <- which.min(substring(pct.L, first = 3))
      vn1 <- pct.L[i]
      disp <- TRUE
    } else vn1 <- "quant"
    if (length(pct.U)) {
      i <- which.min(substring(pct.U, first = 3))
      vn2 <- pct.U[i]
      disp <- TRUE
    } else vn2 <- "quant"
    prec <-  mean(pred[ ,"quant"])/100
    if (disp)  prec <- pmin(prec, min(pred[ ,vn2] - pred[, vn1]), na.rm = TRUE)
    dig.quant <- -floor(log(prec, base = 10))
    if (dig.quant < 0) dig.quant <-  0 
  } 
  ## round all selected variables
  pred.mod <- pred
  ind <- c("quant", pct.L, pct.U)
  pred.mod[ , ind] <- round(pred.mod[ , ind], digits = dig.quant) 
  pred.mod
}

##===================================================================
## add a small amount of noise to x
##
##===================================================================
OTjitter <- function(x, threshold = NULL) {
  d <- diff(sort(x))
  signoise <- pmin(min(d[d>0]) / 5, mean(abs(x))/500)
  mynoise  <- rnorm(length(x), mean = 0, sd = signoise)
  x.noised <- x + mynoise
  if (length(threshold) > 0) {
    if (any(x < threshold)) stop("'threshold' must be a numeric <= min(x)")
    ind <- x.noised < threshold
    x.noised[ind] <- x[ind] + abs(mynoise[ind])
  }
  x.noised  
}
##*****************************************************************************
print.Renouv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L, 
                      quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

##======================================================================
## summary method for 'Renouv'
##======================================================================
summary.Renouv <- function(object,
                           correlation = FALSE,
                           symbolic.cor = FALSE,
                           ...) {
  ## distribution
  ans <- object
  ## improve info on degree of freedom
  ans$df <- as.integer(c(object$df, object$nobs))
  names(ans$df) <- c("par", "obs")
  ## coefficients: take care for fixed ones!
  est <- object$estimate
  se <- sqrt(diag(object$cov))
  tval <- est / se
  ## correction of df due to historical data
  rdf <- object$nobs - object$df
  ans$df <- c("par" = object$df, "obs" = object$nobs, "resid" = rdf)
  ans$coefficients <-
    cbind(est, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <-
    list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$pred <- roundPred(object$pred)
  if (correlation) {
    ans$correlation <- object$corr
    ans$symbolic.cor <- symbolic.cor  
  }
  class(ans) <- "summary.Renouv"
  ans
}

##*****************************************************************************
print.summary.Renouv <-
  function(x,
           coef = TRUE,
           pred = TRUE,
           probT = FALSE,
           digits = max(3, getOption("digits") - 3),
           symbolic.cor = x$symbolic.cor,
           signif.stars = getOption("show.signif.stars"),
           ...) {
  ## cat(sprintf("o Number of OT observations : %d\n", length(x$y.OT)))
  cat(sprintf(paste("o Main sample 'Over Threshold'\n",
                    "   . Threshold        %8.2f\n",
                    "   . Effect. duration %8.2f years\n",
                    "   . Nb. of exceed.   %5d\n\n",
                    collapse = ""),
              x$threshold, x$effDuration, length(x$y.OT)))
  cat(sprintf("o Estimated rate 'lambda' for Poisson process (events): %5.2f evt/year.\n\n",
              x$estimate["lambda"]))
  cat(sprintf("o Distribution for exceedances y: \"%s\", with %d par. ",
              x$distname.y, x$p.y))
  cat(paste(sprintf("\"%s\"", x$parnames.y), collapse = ", "), "\n\n")
  if (x$transFlag){
    cat(sprintf("o Transformation applied: \"%s\"\n\n",
                x$trans.y))
  } else {
    cat("o No transformation applied\n\n")
  }
  if (coef) {
    cat("o Coefficients\n\n")
    if (!probT) {
      print(x$coefficients[ , 1:3])
    } else {
      printCoefmat(x$coefficients,
                   digits = digits,
                   signif.stars = signif.stars,
                   na.print = "NA", ...)
    }
    cat("\n")
    cat(sprintf("Degrees of freedom: %d (param.) and %d (obs)\n", x$df["par"], x$df["obs"]))
    cat("\n")
  }
  if (any(x$fixed)) {
    cat("The following coef. were fixed\n")
    print(names(x$fixed)[x$fixed])
    cat("\n")
  }
  cat(sprintf("o Inference method used for return levels\n\"%s\"\n\n",
              x$infer.method))
  if (pred) {
    cat("o Return levels\n\n")
    print(x$pred)
    cat("\n\n")
  }
  if (x$history.MAX$flag) {
    cat(sprintf("o 'MAX' historical info: %d  blocks, %d obs., total duration = %5.2f years\n\n",
                nlevels(x$history.MAX$block),
                length(unlist(x$history.MAX$data)),
                sum(x$history.MAX$effDuration)))
    if (nlevels(x$history.MAX$block) <= 4L) {
      for (i in 1L:nlevels(x$history.MAX$block) ) {
        Zi <- x$history.MAX$data[[i]]
        ri <- length(Zi)
        if ( ri <= 12L ) {
          obs.str <- paste(format(Zi),  collapse = ", ")
        } else {
          obs.str <- paste(paste(format(Zi[1L:3L]), collapse = ", "),
                           "...",
                           paste(format(Zi[(ri-2L):ri]), collapse = ", "),
                           sep = ", ")
        }
        cat(sprintf("  * block %d, %5.2f years, %d obs.\n\n     %s\n\n",
                    i, x$history.MAX$effDuration[[i]], ri, obs.str))
      }
    }
  } else {
    cat("o no 'MAX' historical data\n\n")
  }
  if (x$history.OTS$flag) {
    cat(sprintf("o 'OTS' historical info: %d  blocks, %d obs., total duration = %5.2f years\n\n",
                nlevels(x$history.OTS$block),
                length(unlist(x$history.OTS$data)),
                sum(x$history.OTS$effDuration)))
    if (nlevels(x$history.OTS$block) <= 4L) {
      for (i in 1L:nlevels(x$history.OTS$block) ) {
        Zi <- x$history.OTS$data[[i]]
        ri <- length(Zi)
        if (ri > 0L) {
          if (ri <= 12L) {
            obs.str <- paste(format(Zi),  collapse = ", ")
          } else {
            obs.str <- paste(paste(format(Zi[1L:3L]), collapse = ", "),
                             "...",
                             paste(format(Zi[(ri-2L):ri]), collapse = ", "),
                             sep = ", ")
          }
          cat(sprintf("  * block %d, %5.2f years, thresh. %8.2f, %d obs.\n\n    %s\n\n",
                      i, x$history.OTS$effDuration[[i]],
                      x$history.OTS$threshold[[i]],
                      ri, obs.str))
          
        } else {
          cat(sprintf("  * block %d, %5.2f years, thresh. %8.2f, %d obs.\n\n",
                      i, x$history.OTS$effDuration[[i]],
                      x$history.OTS$threshold[[i]], ri))
        }
      }
    }
  } else {
    cat("o no 'OTS' historical data\n\n")
  }
  ## 2012-12-22 PREVIOUS CODE
  if (FALSE) {
    if (x$history.MAX$flag) {
      cat(sprintf("o 'MAX' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                  nlevels(x$history.MAX$block),
                  length(unlist(x$history.MAX$data)),
                  sum(x$history.MAX$effDuration))) 
    } else  cat("o no 'MAX' historical data\n\n")
    if (x$history.OTS$flag) {
      cat(sprintf("o 'OTS' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                  nlevels(x$history.OTS$block),
                  length(unlist(x$history.OTS$data)),
                sum(x$history.OTS$effDuration)))
    } else  cat("o no 'OTS' historical data\n\n")
  } 
  cat("o Kolmogorov-Smirnov test\n")
  print(x$KS.test)
  cat("\n")
  ## copied from 'print.summary.lm' in the 'stats' package
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("o Correlation of Coefficients:\n")
      if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
        print(symnum(correl, abbr.colnames = NULL))
        cat("\n")
      } else {
        correl <- format(round(correl, 2), nsmall = 2, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
        cat("\n")
      }
    }
  }
  
  if (!is.null(x$MAX)) {
    cat("o Implied model for block maxima\n")
    cat("  Distribution:", x$MAX$distname, "\n")
    cat("  Coeffficients\n")
    print(x$MAX$estimate)
  }
  
  
}

##===========================================================================
##
##===========================================================================
format.summary.Renouv <-
  function(x,
           ...) {
  ## cat(sprintf("o Number of OT observations : %d\n", length(x$y.OT)))
  text <- sprintf(paste("o Main sample 'Over Threshold'\n",
                    "   . Threshold        %8.2f\n",
                    "   . Effect. duration %8.2f years\n",
                    "   . Nb. of exceed.   %5d\n\n",
                    collapse = ""),
              x$threshold, x$effDuration, length(x$y.OT))
 
  text <- paste(text,
                sprintf("o Estimated rate 'lambda' for Poisson process (events): %5.2f evt/year.\n\n",
                        x$estimate["lambda"]))
 
  text <- paste(text, sprintf("o Distribution for exceedances y: \"%s\", with %d par. ",
                              x$distname.y, x$p.y))

  text <- paste(text, paste(sprintf("\"%s\"", x$parnames.y), collapse = ", "), "\n\n")

  if (any(x$fixed)) {
    cat("The following coef. were fixed\n", names(x$fixed)[x$fixed], "\n")
  }

  if (x$history.MAX$flag) {
    text <- paste(text,
                  sprintf("o 'MAX' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                          nlevels(x$history.MAX$block), length(unlist(x$history.MAX$data)),
                          sum(x$history.MAX$effDuration)))
    
  } else  text <- paste(text, "o no 'MAX' historical data\n\n")
  
  if (x$history.OTS$flag) {
    text <- paste(text, sprintf("o 'OTS' historical info: %d blocks, %d obs., total duration = %5.2f years\n\n",
                                nlevels(x$history.OTS$block), length(unlist(x$history.OTS$data)),
                                sum(x$history.OTS$effDuration)))
  } else  text <- paste(text, "o no 'OTS' historical data\n\n")
  
  text <- paste(text, sprintf("o Kolmogorov-Smirnov test\n   D = %6.4f, p-value = %6.4f\n",
                              x$KS.test$statistic, x$KS.test$p.value))
  text
}

##====================================================================
## coefficients method
##
##====================================================================
coef.Renouv <- function(object, ...) {
  object$estimate
}
##====================================================================
## vcov method
##====================================================================
vcov.Renouv <- function(object, ...) {
  object$cov
}
##====================================================================
## CAUTION.
##
## For now, 'logLik', 'nobs' only works when no historical data are
## used.
##====================================================================
nobs.Renouv <- function(object, ...) {
  object$nobs
}

logLik.Renouv <- function(object, ...) {
  res <- object$logLik
  attr(res, "df") <- object$df
  attr(res, "nobs") <- object$nobs
  res
}
##====================================================================
## AIC methods
##====================================================================
AIC.Renouv <- function(object, ..., k = 2) {
  return(NextMethod())
}

BIC.Renouv <- function(object, ...) {
  return(NextMethod())
}

