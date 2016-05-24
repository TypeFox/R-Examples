ci.reliability <- function(data = NULL, S = NULL, N = NULL, aux = NULL, type = "omega", 
                           interval.type = "ll", B = 10000, conf.level = 0.95) 
{

if(is.null(type) || type == "default") 
    {
    if(!is.null(data) && all(apply(data, 2, function(x) length(table(x))) <= 10)) 
      {type <- 5
      warnings("Categorical omega is used because your variables look like ordered categorical variables. If not, please specify the 'type' argument.")
      } else {
      type <- 4
      warnings("Hiearachical omega is used by default for covariance matrix input or continuous data input.")
              }
  }
  
  type1 <- c(1, "alpha", "true score equivalent", "true-score equivalent", "true score", "equivalent", "tau equivalent", "cronbach", "tau-equivalent", "a")
  type2 <- c(2, "alpha-analytic", "alpha analytic", "alpha factor", "alpha-factor", "alpha cfa", "alpha-cfa", "aa") # Not exactly equal to coefficient alpha because ULS is not used
  type3 <- c(3, "omega", "congeneric", "w")
  type4 <- c(4, "hierarchical omega", "hierarchical", "h")
  type5 <- c(5, "categorical omega", "categorical", "ordered", "ordered categorical", "c", "cat")
  
  type <- tolower(type)
  
  if(type %in% type1) {
    type <- 1
  } else if (type %in% type2) {
    type <- 2
  } else if (type %in% type3) {
    type <- 3
  } else if (type %in% type4) {
    type <- 4
  } else if (type %in% type5) {
    type <- 5
  } else {
    stop("Please provide a correct type of reliability: 'alpha', 'alpha-analytic', 'omega', 'hierarchical', or 'categorical'.")
  }
  
  interval.type <- .translateinterval.type(interval.type, pos = 1)
  
  if (!(is.vector(conf.level) && (length(conf.level) == 1) && is.numeric(conf.level))) 
    stop("Please put a number in the confidence level!")
  
  if(!is.null(data) & !is.null(S)) {
    stop("Both data and covariance matrix cannot be specified simultaneously.")
  } else if (!is.null(S)) {
    return(.ci.reliability.cov(S = S, N = N, type = type, interval.type = interval.type, conf.level = conf.level))
  } else if (!is.null(data)) {
    return(.ci.reliability.data(data = data, aux = aux, type = type, interval.type = interval.type, B = B, conf.level = conf.level))
  } else {
    stop("Either data or covariance matrix must be specified.")
  }
}

.ci.reliability.cov <- function(S = NULL, N = NULL, type = NULL, interval.type = NULL, conf.level = 0.95) {
  if(is.null(N)) stop("Please specify sample size (N)")
  if(!isSymmetric(S, tol = 1e-6)) stop("Covariance matrix must be symmetric.")
  alpha <- 1 - conf.level
  
  varnames <- colnames(S)
  if(is.null(varnames)) {
    varnames <- paste0("y", 1:nrow(S))
    colnames(S) <- rownames(S) <- varnames
  }
  
  ### Point Estimate
  relia <- NA
  se <- NA
  ci.lower <- NA
  ci.upper <- NA
  q <- nrow(S)
  if(type == 1) {
    relia <- .find.alpha(S)
  } else if(type == 2) {
    temp <- .run.cfa.cov(S, N, varnames, se = "none", equal.loading = TRUE)
    relia <- temp$relia
  } else if(type == 3) {
    temp <- .run.cfa.cov(S, N, varnames, se = "none")
    relia <- temp$relia
  } else if(type == 4) {
    relia <- .getH.cov(S, N, varnames)
    if(interval.type != 0) stop("Hierarchical omega requires a full data set to estimate a confidence interval.")
  } else if(type == 5) {
    stop("Categorical omega requires a real data.")
  }
  
  ### Interval Estimate
  crit <- qnorm(1 - (1 - conf.level)/2)
  if(interval.type == 11) {
    if (type == 1) {
      temp <- .interval.type11_type1(relia, q, N, crit)
      se <- temp[1]
      ci.lower <- temp[2]
      ci.upper <- temp[3]
    } else if (type == 2) {
      temp <- .run.cfa.cov(S, N, varnames, se = "default", equal.loading = TRUE, equal.error = TRUE)
      if(!is.na(temp$relia)) relia <- temp$relia # Overwrite because of the parallel assumption
      se <- temp$se
      ci.lower <- relia - (crit * se)
      ci.upper <- relia + (crit * se)
    } else {
      stop("Coefficient omega, hierarchical omega, and categorical omega do not work with the parallel method for interval estimation.")
    }
    
  } else if (interval.type == 12) {
    if (type %in% 1:4) {
      result <- .fMethod(relia, N - 1, (N - 1) * (q - 1), conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else {
      stop("This interval estimation method was not designed for categorical omega.")
    }		
  } else if (interval.type == 13) {
    if (type %in% 1:4) {
      result <- .fMethod(relia, N, N * (q - 1), conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else {
      stop("This interval estimation method was not designed for categorical omega.")
    }		
  } else if (interval.type == 21) {
    result <- .fisherCIrelia(relia, N, crit)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 22) {
    result <- .iccCIrelia(relia, N, q, crit, bonett = TRUE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 23) {
    result <- .hkCIrelia(relia, S, q, N, crit, correct = FALSE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 24) {
    result <- .hkCIrelia(relia, S, q, N, crit, correct = TRUE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 25) {
    result <- .iccCIrelia(relia, N, q, crit, bonett = FALSE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type %in% 31:36) {
    estimator <- "ml"
    logistic <- interval.type %in% c(32, 34, 36)
    if(interval.type %in% 33:34) estimator <- "mlr"
    if(interval.type %in% 35:36) estimator <- "wls"
    if (type == 1) {
      if(interval.type %in% 31:32) {
        result <- .interval.type31_32_type1(relia, S, q, N, crit, logistic = logistic)
        se <- result$se
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
      } else if(interval.type %in% 33:34) {
        stop("MLR estimator is not available for coefficient alpha")
      } else {
        stop("ADF requires a full data set.")		
      }
    } else if (type == 2) {
      result <- .run.cfa.cov(S, N, varnames, se = "default", equal.loading = TRUE)
      if(!is.na(result$relia)) relia <- result$relia # Overwrite it if different estimators are used
      se <- result$se
      if(logistic) {
        temp <- .logisticT(result$relia, result$se, crit)
        ci.upper <- temp[2]
        ci.lower <- temp[1]
      } else {
        ci.lower <- relia - (crit * se)
        ci.upper <- relia + (crit * se)
      }
    } else if (type == 3) {
      result <- .run.cfa.cov(S, N, varnames, se = "default")
      se <- result$se
      if(logistic) {
        temp <- .logisticT(result$relia, result$se, crit)
        ci.upper <- temp[2]
        ci.lower <- temp[1]
      } else {
        ci.lower <- relia - (crit * se)
        ci.upper <- relia + (crit * se)
      }
    } else {
      stop("Any normal-theory or ADF methods are not available for hierarchical omega or categorical omega.")
    }
  } else if (interval.type == 37) {
    if (type == 1) {
      result <- .llCovAlpha(S, N, conf.level = conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else if (type == 2) {
      result <- .llCovOmega(S, N, eqload = TRUE, conf.level = conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else if (type == 3) {
      result <- .llCovOmega(S, N, eqload = FALSE, conf.level = conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper		
    } else {
      stop("The profile-likelihood confidence interval method is not available for hierarchical and categorical omega.")
    }
  } else if (interval.type %in% 41:44) {
    stop("Bootstrap confidence interval requires a real data set.")
  }   	
  
  if (!is.na(ci.lower) && ci.lower < 0) ci.lower <- 0
  if (!is.na(ci.upper) && ci.upper > 1) ci.upper <- 1
  out <- list(est = relia, se = se, ci.lower = ci.lower, ci.upper = ci.upper, conf.level = conf.level, type = c("alpha", "alpha-cfa", "omega", "hierarchical omega", "categorical omega")[type], interval.type = .translateinterval.type(interval.type, pos = 2))
  out
}

.ci.reliability.data <- function(data = NULL, aux = NULL, type = NULL, interval.type = NULL, B = 1000, conf.level = 0.95) {
  if(!is.data.frame(data) & !is.matrix(data)) stop("data must be in a data frame or matrix format.")
  
  if(interval.type %in% c(41, 42, 43, 44, 45) && !(is.vector(B) && (length(B) == 1) && is.numeric(B))) 
    stop("Please put a number in the number of bootstrap argument!")
  
  alpha <- 1 - conf.level
  
  varnames <- colnames(data)
  if(!is.null(aux)) {
    if(is.null(varnames)) stop("'data' should have variable names because the list of auxiliary variables is indicated.")
    if(!all(aux %in% varnames)) stop("Some (or all) auxiliary variables are not in the specified data set.")
    if(type == 5) stop("Categorical CFA does not support the auxiliary variable feature because direct maximum likelihood is not used. Multiple imputation should be used instead.")
    varnames <- setdiff(varnames, aux)
  } else {
    if(is.null(varnames)) {
      varnames <- paste0("y", 1:ncol(data))
      colnames(data) <- varnames
    }
  }
  
  ### Point Estimate
  relia <- NA
  se <- NA
  ci.lower <- NA
  ci.upper <- NA
  N <- nrow(data)
  q <- length(varnames)
  if(type == 1) {
    S <- .findS(data, varnames, aux = aux)
    relia <- .find.alpha(S)
    N <- attr(S, "effn")
  } else if(type == 2) {
    temp <- .run.cfa(data, varnames, aux, estimator = "mlr", missing = "ml", se = "none", equal.loading = TRUE) 
    relia <- temp$relia
    N <- temp$effn
  } else if(type == 3) {
    temp <- .run.cfa(data, varnames, aux, estimator = "mlr", missing = "ml", se = "none")
    relia <- temp$relia
    N <- temp$effn
  } else if(type == 4) {
    relia <- .getH(data, varnames, aux = aux, estimator = "mlr", se = "none", missing = "ml")
    N <- attr(relia, "effn")
    attr(relia, "effn") <- NULL
  } else if(type == 5) {
    relia <- .catOmega(dat = data)
  }
  
  ### Interval Estimate
  crit <- qnorm(1 - (1 - conf.level)/2)
  if(interval.type == 11) {
    if (type == 1) {
      temp <- .interval.type11_type1(relia, q, N, crit)
      se <- temp[1]
      ci.lower <- temp[2]
      ci.upper <- temp[3]
    } else if (type == 2) {
      temp <- .run.cfa(data, varnames, aux, estimator = "mlr", missing = "ml", se = "default", equal.loading = TRUE, equal.error = TRUE)
      if(!is.na(temp$relia)) relia <- temp$relia # Overwrite because of the parallel assumption
      se <- temp$se
      ci.lower <- relia - (crit * se)
      ci.upper <- relia + (crit * se)
    } else {
      stop("Coefficient omega, hierarchical omega, and categorical omega do not work with the parallel method for interval estimation.")
    }
    
  } else if (interval.type == 12) {
    if (type %in% 1:4) {
      result <- .fMethod(relia, N - 1, (N - 1) * (q - 1), conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else {
      stop("This interval estimation method was not designed for categorical omega.")
    }		
  } else if (interval.type == 13) {
    if (type %in% 1:4) {
      result <- .fMethod(relia, N, N * (q - 1), conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else {
      stop("This interval estimation method was not designed for categorical omega.")
    }		
  } else if (interval.type == 21) {
    result <- .fisherCIrelia(relia, N, crit)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 22) {
    result <- .iccCIrelia(relia, N, q, crit, bonett = TRUE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 23) {
    S <- .findS(data, varnames, aux = aux)
    result <- .hkCIrelia(relia, S, q, N, crit, correct = FALSE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 24) {
    S <- .findS(data, varnames, aux = aux)
    result <- .hkCIrelia(relia, S, q, N, crit, correct = TRUE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type == 25) {
    result <- .iccCIrelia(relia, N, q, crit, bonett = FALSE)
    se <- result$se
    ci.lower <- result$ci.lower
    ci.upper <- result$ci.upper
  } else if (interval.type %in% 31:36) {
    estimator <- "ml"
    logistic <- interval.type %in% c(32, 34, 36)
    if(interval.type %in% 33:34) estimator <- "mlr"
    if(interval.type %in% 35:36) estimator <- "wls"
    if (type == 1) {
      if(interval.type %in% 31:32) {
        S <- .findS(data, varnames, aux = aux)
        result <- .interval.type31_32_type1(relia, S, q, N, crit, logistic = logistic)
        se <- result$se
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
      } else if(interval.type %in% 33:34) {
        stop("MLR estimator is not available for coefficient alpha")
      } else {
        if(any(is.na(data))) stop("data should not have any missing values to use this adf method.")
        se <- .seReliabilityAdf(data)
        if(logistic) {
          temp <- .logisticT(relia, se, crit)
          ci.upper <- temp[2]
          ci.lower <- temp[1]
        } else {
          ci.lower <- relia - (crit * se)
          ci.upper <- relia + (crit * se)
        }			
      }
    } else if (type == 2) {
      result <- .run.cfa(data, varnames, aux, estimator = estimator, missing = "ml", se = "default", equal.loading = TRUE)
      if(!is.na(result$relia)) relia <- result$relia # Overwrite it if different estimators are used
      se <- result$se
      if(logistic) {
        temp <- .logisticT(result$relia, result$se, crit)
        ci.upper <- temp[2]
        ci.lower <- temp[1]
      } else {
        ci.lower <- relia - (crit * se)
        ci.upper <- relia + (crit * se)
      }
    } else if (type == 3) {
      result <- .run.cfa(data, varnames, aux, estimator = estimator, missing = "ml", se = "default")
      se <- result$se
      if(logistic) {
        temp <- .logisticT(result$relia, result$se, crit)
        ci.upper <- temp[2]
        ci.lower <- temp[1]
      } else {
        ci.lower <- relia - (crit * se)
        ci.upper <- relia + (crit * se)
      }
    } else {
      stop("Any normal-theory or ADF methods are not available for hierarchical omega or categorical omega.")
    }
  } else if (interval.type == 37) {
    if(!is.null(aux)) stop("The profile-likelihood method does not support the auxiliary-variable feature.")
    if (type == 1) {
      result <- .llDataAlpha(data, conf.level = conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else if (type == 2) {
      result <- .llDataOmega(data, eqload = TRUE, conf.level = conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper
    } else if (type == 3) {
      result <- .llDataOmega(data, eqload = FALSE, conf.level = conf.level)
      ci.lower <- result$ci.lower
      ci.upper <- result$ci.upper		
    } else {
      stop("The profile-likelihood confidence interval method is not available for hierarchical and categorical omega.")
    }
  } else if (interval.type %in% 41:44) {
    boot.out <- NULL
    .bs1 <- function(data, i, varnames) {
      S <- cov(data[i, varnames])
      .find.alpha(S)
    }
    .bs1miss <- function(data, i, varnames, aux) {
      S <- .findS(data[i,], varnames = varnames, aux = aux)
      .find.alpha(S)
    }
    .bs2 <- function(data, i, varnames, aux) {
      temp <- .run.cfa(data[i,], varnames = varnames, aux = aux, estimator = "mlr", missing = "ml", se = "none", equal.loading = TRUE) 
      relia <- temp$relia
    }
    .bs3 <- function(data, i, varnames, aux) {
      temp <- .run.cfa(data[i,], varnames = varnames, aux = aux, estimator = "mlr", missing = "ml", se = "none", equal.loading = FALSE) 
      relia <- temp$relia
    }
    .bs4 <- function(data, i, varnames, aux) {
      result <- .getH(data[i,], varnames = varnames, aux = aux, estimator = "mlr", se = "none", missing = "ml")
      attr(result, "effn") <- NULL
      result
    }
    .bs5 <- function(data, i, varnames) {
      .catOmega(dat = data[i, varnames])
    }	
    if (type == 1) {
      if(any(is.na(data))) {
        boot.out <- boot::boot(data = data, statistic = .bs1miss, R = B, stype = "i", varnames = varnames, aux = aux)
      } else {
        boot.out <- boot::boot(data = data, statistic = .bs1, R = B, stype = "i", varnames = varnames)
      }
    } else if (type == 2) {
      boot.out <- boot::boot(data = data, statistic = .bs2, R = B, stype = "i", varnames = varnames, aux = aux)
    } else if (type == 3) {
      boot.out <- boot::boot(data = data, statistic = .bs3, R = B, stype = "i", varnames = varnames, aux = aux)
    } else if (type == 4) {
      boot.out <- boot::boot(data = data, statistic = .bs4, R = B, stype = "i", varnames = varnames, aux = aux)
    } else if (type == 5) {
      boot.out <- boot::boot(data = data, statistic = .bs5, R = B, stype = "i", varnames = varnames)
    } else {
      stop("Something was wrong at the if..else in bootstrap.")
    }
    ci.output <- NULL
    se <- apply(boot.out$t, 2, sd, na.rm = TRUE)
    if(interval.type == 41) {
      ci.lower <- relia - (crit * se)
      ci.upper <- relia + (crit * se)
    } else if (interval.type == 42) {
      temp <- .logisticT(relia, se, crit)
      ci.upper <- temp[2]
      ci.lower <- temp[1]
    } else if (interval.type == 43) {
      ci.output <- boot::boot.ci(boot.out = boot.out, conf = conf.level, type = "perc")$perc
      ci.lower <- ci.output[4]
      ci.upper <- ci.output[5]
    } else if (interval.type == 44) {
      if (B < nrow(data)) 
        warnings("The number of bootstrap samples is less than the number of cases.\nIf 'bca' cannot be calculated, please make sure that the number of cases is greater than\nthe number of bootstrap samples.")
      # From https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html
      ci.output <- boot::boot.ci(boot.out = boot.out, conf = conf.level, type = "bca")$bca
      ci.lower <- ci.output[4]
      ci.upper <- ci.output[5]
    }
  }   	
  if (!is.na(ci.lower) && ci.lower < 0) ci.lower <- 0
  if (!is.na(ci.upper) && ci.upper > 1) ci.upper <- 1
  out <- list(est = relia, se = se, ci.lower = ci.lower, ci.upper = ci.upper, conf.level = conf.level, type = c("alpha", "alpha-cfa", "omega", "hierarchical omega", "categorical omega")[type], interval.type = .translateinterval.type(interval.type, pos = 2))
  out
}

.find.alpha <- function(S = NULL) {
  q <- nrow(S)
  sigma.jj <- sum(diag(S))
  sigma2.Y <- sum(S)
  (q/(q - 1)) * (1 - sigma.jj/sigma2.Y)
}

.findS <- function(data, varnames, aux = NULL) {
  script <- NULL
  for(i in 2:length(varnames)) {
    script <- c(script, paste(varnames[i], "~~", paste("NA*", varnames[1:i], collapse = "+")))
  }
  fit <- NULL
  if(!is.null(aux)) {
    try(fit <- semTools::cfa.auxiliary(script, data = data, aux = aux, std.lv = TRUE, missing = "ml", se = "none"), silent = TRUE)
  } else {
    try(fit <- lavaan::cfa(script, data = data, std.lv = TRUE, missing = "ml", se = "none"), silent = TRUE)
  }
  S <- lavaan::inspect(fit, "cov.ov")[varnames, varnames]
  pe <- lavaan::parameterEstimates(fit)
  N <- nrow(data)
  if("fmi" %in% colnames(pe)) {
    # Effective Sample Size = N(1 - \lambda) where \lambda = average FMI
    fmi <- pe[,"fmi"]
    N <- N * (1 - mean(fmi, na.rm = TRUE))
  }
  attr(S, "effn") <- ceiling(N)
  S
}

.run.cfa <- function(data, varnames, aux, estimator = "default", se = "default", missing = "default", equal.loading = FALSE, equal.error = FALSE) {
  q <- length(varnames)
  N <- nrow(data)
  if (equal.loading) {
    loadingName <- rep("a1", q)
  } else {
    loadingName <- paste("a", 1:q, sep = "")
  }
  if (equal.error) {
    errorName <- rep("b1", q)
  } else {
    errorName <- paste("b", 1:q, sep = "")
  }
  model <- paste0("f1 =~ NA*", varnames[1], " + ")
  loadingLine <- paste(paste(loadingName, "*", varnames, sep = ""), collapse = " + ")
  factorLine <- "f1 ~~ 1*f1\n"
  errorLine <- paste(paste(varnames, " ~~ ", errorName, "*", varnames, sep = ""), collapse = "\n")
  sumLoading <- paste("loading :=", paste(loadingName, collapse = " + "), "\n")
  sumError <- paste("error :=", paste(errorName, collapse = " + "), "\n")
  relia <- "relia := (loading^2) / ((loading^2) + error) \n"
  model <- paste(model, loadingLine, "\n", factorLine, errorLine, "\n", sumLoading, sumError, relia)
  if(!is.null(aux)) {
    e <- try(fit <- semTools::cfa.auxiliary(model, data = data, aux = aux, missing = missing, se = se, estimator = estimator), silent = TRUE)
  } else {
    e <- try(fit <- lavaan::cfa(model, data = data, missing = missing, se = se, estimator = estimator), silent = TRUE)
  }
  converged <- FALSE
  if(is(e, "try-error")) {
    converged <- FALSE
  } else {
    converged <- fit@Fit@converged
    errorcheck <- diag(lavaan::inspect(fit, "se")$theta)
    if (se != "none" && any(errorcheck <= 0)) converged <- FALSE		
  }
  if (converged) {
    loading <- unique(as.vector(lavaan::inspect(fit, "coef")$lambda))
    error <- unique(diag(lavaan::inspect(fit, "se")$theta))
    pe <- lavaan::parameterEstimates(fit)
    r <- which(pe[,"lhs"] == "relia")
    u <- pe[which(pe[,"lhs"] == "loading"), "est"]
    v <- pe[which(pe[,"lhs"] == "error"), "est"]
    est <- pe[r, "est"]
    if (se == "none") {
      paramCov <- NULL
      stderr <- NA
    } else {
      paramCov <- lavaan::vcov(fit)
      stderr <- pe[r, "se"]
    }
    if("fmi" %in% colnames(pe)) {
      # Effective Sample Size = N(1 - \lambda) where \lambda = average FMI
      fmi <- pe[,"fmi"]
      N <- N * (1 - mean(fmi, na.rm = TRUE))
    }
  } else {
    loading <- NA
    error <- NA
    if (se == "none") {
      paramCov <- NULL
    } else {
      paramCov <- NA
    }
    u <- NA
    v <- NA
    est <- NA
    stderr <- NA
  }
  result <- list(load = loading, error = error, 
                 vcov = paramCov, converged = converged, u = u, v = v, relia = est, se = stderr, effn = ceiling(N))
  return(result)
}

.run.cfa.cov <- function(S, N, varnames, se = "default", equal.loading = FALSE, equal.error = FALSE) {
  q <- length(varnames)
  if (equal.loading) {
    loadingName <- rep("a1", q)
  } else {
    loadingName <- paste("a", 1:q, sep = "")
  }
  if (equal.error) {
    errorName <- rep("b1", q)
  } else {
    errorName <- paste("b", 1:q, sep = "")
  }
  model <- paste0("f1 =~ NA*", varnames[1], " + ")
  loadingLine <- paste(paste(loadingName, "*", varnames, sep = ""), collapse = " + ")
  factorLine <- "f1 ~~ 1*f1\n"
  errorLine <- paste(paste(varnames, " ~~ ", errorName, "*", varnames, sep = ""), collapse = "\n")
  sumLoading <- paste("loading :=", paste(loadingName, collapse = " + "), "\n")
  sumError <- paste("error :=", paste(errorName, collapse = " + "), "\n")
  relia <- "relia := (loading^2) / ((loading^2) + error) \n"
  model <- paste(model, loadingLine, "\n", factorLine, errorLine, "\n", sumLoading, sumError, relia)
  e <- try(fit <- lavaan::cfa(model, sample.cov = S, estimator = "default", sample.nobs = N, se = se), silent = TRUE)
  
  converged <- FALSE
  if(is(e, "try-error")) {
    converged <- TRUE
  } else {
    converged <- fit@Fit@converged
    errorcheck <- diag(lavaan::inspect(fit, "se")$theta)
    if (se != "none" && any(errorcheck <= 0)) converged <- FALSE		
  }
  if (converged) {
    loading <- unique(as.vector(lavaan::inspect(fit, "coef")$lambda))
    error <- unique(diag(lavaan::inspect(fit, "se")$theta))
    pe <- lavaan::parameterEstimates(fit)
    r <- which(pe[,"lhs"] == "relia")
    u <- pe[which(pe[,"lhs"] == "loading"), "est"]
    v <- pe[which(pe[,"lhs"] == "error"), "est"]
    est <- pe[r, "est"]
    if (se == "none") {
      paramCov <- NULL
      stderr <- NA
    } else {
      paramCov <- lavaan::vcov(fit)
      stderr <- pe[r, "se"]
    }
    if("fmi" %in% colnames(pe)) {
      # Effective Sample Size = N(1 - \lambda) where \lambda = average FMI
      fmi <- pe[,"fmi"]
      N <- N * (1 - mean(fmi, na.rm = TRUE))
    }
  } else {
    loading <- NA
    error <- NA
    if (se == "none") {
      paramCov <- NULL
    } else {
      paramCov <- NA
    }
    u <- NA
    v <- NA
    est <- NA
    stderr <- NA
  }
  result <- list(load = loading, error = error, 
                 vcov = paramCov, converged = converged, u = u, v = v, relia = est, se = stderr, effn = ceiling(N))
  return(result)
}

.getH <- function(data, varnames, aux, estimator = "default", se = "default", missing = "default") {
  result <- .run.cfa(data = data, varnames = varnames, aux = aux, estimator = estimator, se = se, missing = missing, equal.loading = FALSE, equal.error = FALSE)
  S <- .findS(data, varnames, aux = aux)
  h <- (result$u)^2/sum(S)
  attr(h, "effn") <- result$effn
  h
}

.getH.cov <- function(S, N, varnames) {	
  result <- .run.cfa.cov(S = S, N = N, varnames = varnames, se = "none")
  h <- (result$u)^2/sum(S)
  h
}

.catOmega <- function(dat) {
  
if(!requireNamespace("lavaan", quietly = TRUE)) stop("The package 'lavaan' is needed; please install the package and try again.")
if(!requireNamespace("mnormt", quietly = TRUE)) stop("The package 'mnormt' is needed; please install the package and try again.")

  q <- ncol(dat)
  for(i in 1:q) dat[,i] <- ordered(dat[,i])
  varnames <- paste0("y", 1:q)
  colnames(dat) <- varnames
  loadingName <- paste("a", 1:q, sep = "")
  errorName <- paste("b", 1:q, sep = "")
  model <- "f1 =~ NA*y1 + "
  loadingLine <- paste(paste(loadingName, "*", varnames, sep = ""), collapse = " + ")
  factorLine <- "f1 ~~ 1*f1\n"
  model <- paste(model, loadingLine, "\n", factorLine)
  error <- try(fit <- lavaan::cfa(model, data = dat, se = "none", ordered = varnames), silent = TRUE)
  converged <- FALSE
  if(!is(error, "try-error") && fit@Fit@converged) converged <- TRUE
  reliab <- NA
  if(converged) {
    param <- lavaan::inspect(fit, "coef")
    ly <- param[["lambda"]]
    ps <- param[["psi"]]
    truevar <- ly%*%ps%*%t(ly)
    threshold <- .getThreshold(fit)[[1]]
    denom <- .polycorLavaan(fit, dat)[varnames, varnames]
    invstdvar <- 1 / sqrt(diag(fit@Fit@Sigma.hat[[1]]))
    polyr <- diag(invstdvar) %*% truevar %*% diag(invstdvar)
    sumnum <- 0
    addden <- 0
    for(j in 1:q) {
      for(jp in 1:q) {
        sumprobn2 <- 0
        addprobn2 <- 0
        t1 <- threshold[[j]]
        t2 <- threshold[[jp]]
        for(c in 1:length(t1)) {
          for(cp in 1:length(t2)) {
            sumprobn2 <- sumprobn2 + .p2(t1[c], t2[cp], polyr[j, jp])
            addprobn2 <- addprobn2 + .p2(t1[c], t2[cp], denom[j, jp])
          }
        }
        sumprobn1 <- sum(pnorm(t1))
        sumprobn1p <- sum(pnorm(t2))
        sumnum <- sumnum + (sumprobn2 - sumprobn1 * sumprobn1p)
        addden <- addden + (addprobn2 - sumprobn1 * sumprobn1p)
      }
    }
    reliab <- sumnum / addden
  }
  reliab
}

.p2 <- function(t1, t2, r) {
  mnormt::pmnorm(c(t1, t2), c(0,0), matrix(c(1, r, r, 1), 2, 2))
}


.polycorLavaan <- function(object, data) {
  ngroups <- object@Data@ngroups
  coef <- lavaan::inspect(object, "coef")
  targettaunames <- NULL
  if(ngroups == 1) {
    targettaunames <- rownames(coef$tau)
  } else {
    targettaunames <- rownames(coef[[1]]$tau)
  }
  barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
  varnames <- unique(apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2])))
  script <- ""
  for(i in 2:length(varnames)) {
    temp <- paste0(varnames[1:(i - 1)], collapse = " + ")
    temp <- paste0(varnames[i], "~~", temp, "\n")
    script <- paste(script, temp)
  }
  suppressWarnings(newobject <- .refit(script, data, varnames, object))
  if(ngroups == 1) {
    return(lavaan::inspect(newobject, "coef")$theta)
  } else {
    return(lapply(lavaan::inspect(newobject, "coef"), "[[", "theta"))
  }
}

.getThreshold <- function(object) {
  ngroups <- object@Data@ngroups
  coef <- lavaan::inspect(object, "coef")
  result <- NULL
  if(ngroups == 1) {
    targettaunames <- rownames(coef$tau)
    barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
    varthres <- apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2]))
    result <- list(split(coef$tau, varthres))
  } else {
    result <- list()
    for(g in 1:ngroups) {
      targettaunames <- rownames(coef[[g]]$tau)
      barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
      varthres <- apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2]))
      result[[g]] <- split(coef[[g]]$tau, varthres)
    }
  }
  return(result)
}


.refit <- function(pt, data, vnames, object) {
  previousCall <- object@call
  args <- as.list(previousCall[-1])
  args$model <- pt
  args$data <- data
  args$ordered <- vnames
  funcall <- as.character(previousCall[[1]])
  tempfit <- do.call(funcall[length(funcall)], args)
}

.interval.type11_type1 <- function(relia, q, N, crit) {
  variance <- (2 * (1 - relia)^2 * q)/((N - 1) * (q - 1))
  se <- as.vector(sqrt((2 * (1 - relia)^2 * q)/((N - 1) * (q - 1))))
  c(se, relia - (crit * se), relia + (crit * se))
}

.interval.type31_32_type1 <- function(relia, S, q, N, crit, logistic = FALSE) {
  cor.mat <- cov2cor(S)
  j <- cbind(rep(1, times = q))
  step.1 <- (q^2/(q - 1)^2)
  gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
  gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
  gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
  gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
  gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% cor.mat) %*% j)
  gamma.2 <- gamma.2.1 - gamma.2.2
  gamma.final <- gamma.1 * gamma.2
  variance <- (step.1 * gamma.final)/(N - 1)
  se <- as.vector(sqrt(variance))
  if(logistic) {
    temp <- .logisticT(relia, se, crit)
    ci.upper <- temp[2]
    ci.lower <- temp[1]
  } else {
    ci.upper <- relia + (crit * se)
    ci.lower <- relia - (crit * se)
  }
  return(list(se = se, ci.lower = ci.lower, ci.upper = ci.upper))
}

.fMethod <- function(relia, df1, df2, conf.level = 0.95) {
  alpha <- 1 - conf.level
  fa <- qf(alpha/2, df1, df2)
  fb <- qf(1 - alpha/2, df1, df2)
  ci.lower <- 1 - ((1 - relia) * fb)
  ci.upper <- 1 - ((1 - relia) * fa)
  return(list(ci.lower = ci.lower, ci.upper = ci.upper))
}

.fisherCIrelia <- function(relia, n, crit) {
  z <- 0.5 * log((1 + relia)/(1 - relia))
  se <- sqrt(1/(n - 3))
  me <- se * crit
  zlower <- z - me
  zupper <- z + me
  ci.lower <- (exp(2 * zlower) - 1)/(exp(2 * zlower) + 1)
  ci.upper <- (exp(2 * zupper) - 1)/(exp(2 * zupper) + 1)
  try(if (ci.lower < 0) 
    ci.lower = 0)
  try(if (ci.upper > 1) 
    ci.upper = 1)
  return(list(z = z, se = se, ci.lower = ci.lower, ci.upper = ci.upper))
}

.iccCIrelia <- function(relia, n, q, crit, bonett = TRUE) {
  est <- log(1 - relia)
  se <- NULL
  if (bonett) {
    se <- sqrt(2 * q/((q - 1) * (n - 2)))
  } else {
    se <- sqrt(2 * q/((q - 1) * n))
  }
  me <- se * crit
  ci.lower <- 1 - exp(est + me)
  ci.upper <- 1 - exp(est - me)
  try(if (ci.lower < 0) 
    ci.lower = 0)
  try(if (ci.upper > 1) 
    ci.upper = 1)
  return(list(z = est, se = se, ci.lower = ci.lower, ci.upper = ci.upper))
}

.hkCIrelia <- function(relia, S, q, N, crit, correct = FALSE) {
  m <- 1
  if (correct) {
    averageCov <- mean(S[lower.tri(S, diag = FALSE)])
    Snew <- matrix(averageCov, q, q)
    diag(Snew) <- mean(diag(S))
    lstar <- det(S)/det(Snew)
    l <- (-N * log(lstar))/((q^2 + q - 4)/2)
    cstar <- 1.452 - (0.464 * l) + (0.046 * l^2)
    m <- min(cstar, 1)
  }
  dfn <- (N - 1) * m
  dfd <- (N - 1) * (q - 1) * m
  a <- (1 - (2/(9 * dfn)))/(1 - (2/(9 * dfd)))
  sigma.num <- (2 * q/(9 * dfd)) * ((1 - relia)^(2/3))
  sigma.denom <- (1 - (2/(9 * dfn)))^2
  se <- sqrt(sigma.num/sigma.denom)
  transform.relia <- (1 - relia)^(1/3)
  transform.lower <- transform.relia + (crit * se)
  transform.upper <- transform.relia - (crit * se)
  ci.lower <- 1 - ((a^3) * (transform.lower^3))
  ci.upper <- 1 - ((a^3) * (transform.upper^3))
  return(list(se = se, ci.lower = ci.lower, ci.upper = ci.upper))
}

.logisticT <- function(est, se, crit) {
  if(is.na(se)) return(c(NA, NA))
  logest <- log(est/(1 - est))
  logse <- se / (est * (1 - est))
  loglower <- logest - crit * logse
  logupper <- logest + crit * logse
  if(logupper < loglower) {
    temp <- loglower
    loglower <- logupper
    loguppper <- temp
  }
  lower <- 1 / (1 + exp(-loglower))
  upper <- 1 / (1 + exp(-logupper))
  c(lower, upper)
}

.seReliabilityAdf <- function(data) {
  data <- na.omit(data)
  N <- dim(data)[1]
  q <- dim(data)[2]
  S <- cov(data)
  vecs.S <- NULL
  for (i in 1:q) {
    vecs.S <- c(vecs.S, S[i:q, i])
  }
  item.mean <- apply(data, 2, mean)
  off.diag.1 <- ((2 * q)/(q - 1))
  off.diag.2 <- (sum(diag(S))/(sum(S)^2))
  off.diag <- off.diag.1 * off.diag.2
  delta <- matrix(rep(off.diag, (q * q)), nrow = q)
  on.diag.1 <- ((-q)/(q - 1))
  on.diag.2 <- ((sum(S) - sum(diag(S)))/(sum(S)^2))
  on.diag <- on.diag.1 * on.diag.2
  diag(delta) <- on.diag
  est.delta <- NULL
  for (i in 1:q) {
    est.delta <- c(est.delta, delta[i:q, i])
  }
  add <- 0
  data <- scale(data, center = TRUE, scale = FALSE)
  for (i in 1:N) {
    dev.y <- as.matrix(data[i,])
    temp <- dev.y %*% t(dev.y)
    si <- NULL
    for (j in 1:q) {
      si <- c(si, temp[j:q, j])
    }
    d.si <- si - vecs.S
    product <- t(est.delta) %*% d.si %*% t(d.si) %*% est.delta
    add <- add + product
  }
  result <- add/(N * (N - 1))
  return(as.vector(sqrt(result)))
}


# Must have for Type 1, 2, and 3 -- S and data
.llCovAlpha <- function(S, N, conf.level = 0.95) {
  q <- dim(S)[1]
  vecs.S <- NULL
  for (i in 1:q) {
    vecs.S <- c(vecs.S, S[1:q, i])
  }
  varnames <- colnames(S)
  if (is.null(varnames)) {
    for (i in 1:q) {
      temp <- paste("x", i, sep = "")
      varnames <- c(varnames, temp)
    }
  }
  # finding Alpha
  C <- NULL
  matrixC <- OpenMx::mxMatrix(type = "Symm", nrow = q, ncol = q, free = TRUE, values = vecs.S, 
                              name = "C")
  matrixP <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = q, 
                              name = "q")
  ONE <- OpenMx::mxMatrix(type = "Full", nrow = q, ncol = 1, free = FALSE, values = 1, name = "ONE")
  ALPHA <- OpenMx::mxAlgebra(expression = (q/(q - 1)) * (1 - (OpenMx::tr(C)/sum(C))), name = "Alpha")
  ALPHAObj <- OpenMx::mxExpectationNormal(covariance = "C", dimnames = varnames)
  Data <- OpenMx::mxData(observed = S, type = "cov", numObs = N)
  ALPHAci <- OpenMx::mxCI("Alpha", interval = conf.level)
  AlphaModel <- OpenMx::mxModel("FindAlpha", matrixC, matrixP, ALPHA, ALPHAObj, Data, ALPHAci, OpenMx::mxFitFunctionML(), ONE)
  AlphaFit <- OpenMx::mxRun(AlphaModel, intervals = T)
  result <- AlphaFit@output$confidenceIntervals[c(1,3)]
  return(list(ci.lower = result[1], ci.upper = result[2]))
} 

# Must have for Type 1, 2, and 3 -- S and data
.llCovOmega <- function(S, N, eqload = FALSE, conf.level = 0.95) {
  q <- dim(S)[1]
  vecs.S <- NULL
  for (i in 1:q) {
    vecs.S <- c(vecs.S, S[1:q, i])
  }
  varnames <- colnames(S)
  if (is.null(varnames)) {
    for (i in 1:q) {
      temp <- paste("x", i, sep = "")
      varnames <- c(varnames, temp)
    }
  }
  Data <- OpenMx::mxData(observed = S, type = "cov", numObs = N)
  Alab <- NA
  A <- U <- L <- NULL # Added to appease CRAN via check (these are only OpenMx variables; calling on them within mxMatrix (without this line) makes the check say "no visible binding for global variable." 
  if(eqload) Alab <- "loadeq"
  matrixA <- OpenMx::mxMatrix(type = "Full", nrow = q, ncol = 1, free = TRUE, values = 0.5, labels = Alab, name = "A")
  matrixL <- OpenMx::mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = 1, 
                              name = "L")
  matrixU <- OpenMx::mxMatrix(type = "Diag", nrow = q, ncol = q, free = TRUE, values = 0.2, 
                              lbound = 1e-04, name = "U")
  OMEGA <- OpenMx::mxAlgebra(expression = (sum(A) * sum(A))/((sum(A) * sum(A)) + sum(U)), 
                             name = "Omega")
  algebraR <- OpenMx::mxAlgebra(expression = A %*% L %*% t(A) + U, name = "R")
  OMEGAObj <- OpenMx::mxExpectationNormal(covariance = "R", dimnames = varnames)
  OMEGAci <- OpenMx::mxCI("Omega", interval = conf.level)
  OmegaModel <- OpenMx::mxModel("FindOmega", matrixA, matrixL, matrixU, OMEGA, algebraR, 
                                OMEGAObj, Data, OMEGAci, OpenMx::mxFitFunctionML())
  OmegaFit <- OpenMx::mxRun(OmegaModel, intervals = T)
  result <- OmegaFit@output$confidenceIntervals[c(1,3)]
  return(list(ci.lower = result[1], ci.upper = result[2]))
} 

.llDataAlpha <- function(data, conf.level = 0.95) {
  q <- ncol(data)
  varnames <- colnames(data)
  S <- cov(data, use = "pairwise.complete.obs")
  vecs.S <- NULL
  for (i in 1:q) {
    vecs.S <- c(vecs.S, S[1:q, i])
  }
  Data <- OpenMx::mxData(observed = data, type = "raw")
  
  # finding Alpha
  C <- NULL
  matrixC <- OpenMx::mxMatrix(type = "Symm", nrow = q, ncol = q, free = TRUE, values = vecs.S, name = "C")
  matrixP <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = q, name = "q")
  ALPHA <- OpenMx::mxAlgebra(expression = (q/(q - 1)) * (1 - (OpenMx::tr(C)/sum(C))), name = "Alpha")
  MEAN <- OpenMx::mxMatrix(type="Full", nrow=1, ncol = q, values = colMeans(data, na.rm = TRUE), free = TRUE, name="M")
  ALPHAObj <- OpenMx::mxExpectationNormal(covariance = "C", means = "M", dimnames = varnames)
  ALPHAci <- OpenMx::mxCI("Alpha", interval = conf.level)
  AlphaModel <- OpenMx::mxModel("FindAlpha", matrixC, matrixP, ALPHA, ALPHAObj, Data, ALPHAci, MEAN,  OpenMx::mxFitFunctionML())
  AlphaFit <- OpenMx::mxRun(AlphaModel, intervals = T)
  result <- AlphaFit@output$confidenceIntervals[c(1,3)]
  return(list(ci.lower = result[1], ci.upper = result[2]))
} 

.llDataOmega <- function(data, eqload = FALSE, conf.level = 0.95) {
  q <- ncol(data)
  varnames <- colnames(data)
  S <- cov(data, use = "pairwise.complete.obs")
  vecs.S <- NULL
  for (i in 1:q) {
    vecs.S <- c(vecs.S, S[1:q, i])
  }
  Data <- OpenMx::mxData(observed = data, type = "raw")
  
  # finding Alpha
  A <- U <- L <- NULL # Added to appease CRAN via check (these are only OpenMx variables; calling on them within mxMatrix (without this line) makes the check say "no visible binding for global variable." 
  MEAN <- OpenMx::mxMatrix(type="Full", nrow=1, ncol = q, values = colMeans(data, na.rm = TRUE), free = TRUE, name="M")
  Alab <- NA
  if(eqload) Alab <- "loadeq"
  # finding Omega
  matrixA <- OpenMx::mxMatrix(type = "Full", nrow = q, ncol = 1, free = TRUE, values = 0.5, labels = Alab, name = "A")
  matrixL <- OpenMx::mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = 1, name = "L")
  matrixU <- OpenMx::mxMatrix(type = "Diag", nrow = q, ncol = q, free = TRUE, values = 0.2, lbound = 1e-04, name = "U")
  OMEGA <- OpenMx::mxAlgebra(expression = (sum(A) * sum(A))/((sum(A) * sum(A)) + sum(U)), 
                             name = "Omega")
  algebraR <- OpenMx::mxAlgebra(expression = A %*% L %*% t(A) + U, name = "R")
  OMEGAObj <- OpenMx::mxExpectationNormal(covariance = "R", mean = "M", dimnames = varnames)
  OMEGAci <- OpenMx::mxCI("Omega", interval = conf.level)
  OmegaModel <- OpenMx::mxModel("FindOmega", matrixA, matrixL, matrixU, OMEGA, algebraR, 
                                OMEGAObj, Data, OMEGAci, MEAN, OpenMx::mxFitFunctionML())
  OmegaFit <- OpenMx::mxRun(OmegaModel, intervals = T)
  result <- OmegaFit@output$confidenceIntervals[c(1,3)]
  return(list(ci.lower = result[1], ci.upper = result[2]))
} 

.translateinterval.type <- function(interval.type, pos = 1) {
  
  # None
  interval.type0 <- c(0, "none", "na")
  
  # Analytic
  interval.type11 <- c(11, "parallel", "sb") # Type 1
  interval.type12 <- c(12, "feldt", "feldt65", "feldt1965", "f", "fdist", 1) # Type 9
  interval.type13 <- c(13, "siotani", "shf", "shf85", "siotani85", "siotani1985", "f2") # Type 12
  
  # Transformation
  interval.type21 <- c(21, "fisher", "naivefisher") # Type 4
  interval.type22 <- c(22, "bonett", 2) # Type 5
  interval.type23 <- c(23, "hakstian-whalen", "hakstianwhalen", "hw", "hakstian", "hw76", "hw1976", "cuberoot") # Type 10
  interval.type24 <- c(24, "hakstian-barchard", "hakstianbarchard", "hb", "hb00", "hb2000", "randomitem") # Type 11
  interval.type25 <- c(25, "intraclass correlation", "icc", "modifiedfisher", "improvedfisher") # Type 13
  
  # ML
  interval.type31 <- c(31, "maximum likelihood (wald ci)", "ml", "normal-theory", "wald", "normal") # Type 2
  interval.type32 <- c(32, "maximum likelihood (logistic ci)", "mll", "logistic", "normall", "waldl", "normal-theory-l") # Type 2 + logistic
  interval.type33 <- c(33, "robust maximum likelihood (wald ci)", "mlr", "robust", "robust ml") # Type 2 + mlr
  interval.type34 <- c(34, "robust maximum likelihood (logistic ci)", "mlrl", "robustl", "robust mll", 3) # Type 2 + logistic + mlr
  interval.type35 <- c(35, "asymptotic distribution free (wald ci)", "adf", "wls") # Type 3
  interval.type36 <- c(36, "asymptotic distribution free (logistic ci)", "adfl", "wlsl") # Type 3 + logistic
  interval.type37 <- c(37, "profile-likelihood", "ll", "likelihood", "logl") # Type 6
  
  # Bootstrap
  interval.type41 <- c(41, "bootstrap standard error (wald ci)", "bsi", "bse", "standardboot", "sdboot") # Type 14
  interval.type42 <- c(42, "bootstrap standard error (logistic ci)", "bsil", "bsel", "standardbootl", "sdbootl") # Type 14 + logistic
  interval.type43 <- c(43, "percentile bootstrap", "perc", "percentile", "percentile ci") # Type 7
  interval.type44 <- c(44, "bca bootstrap", "bca", "boot", "bootstrap", "bias-corrected and acceleration", 4) # Type 8
  
  if(interval.type == "default") {
    interval.type <- 0
    warnings("No confidence interval method is specified. Please specify the confidence interval method in the 'interval.type' argument if you wish to have one. For example, interval.type = 'boot'")
  }
  
  if(interval.type %in% interval.type0) {
    interval.type <- interval.type0[pos]
  } else if (interval.type %in% interval.type11) {
    interval.type <- interval.type11[pos]
  } else if (interval.type %in% interval.type12) {
    interval.type <- interval.type12[pos]
  } else if (interval.type %in% interval.type13) {
    interval.type <- interval.type13[pos]
  } else if (interval.type %in% interval.type21) {
    interval.type <- interval.type21[pos]
  } else if (interval.type %in% interval.type22) {
    interval.type <- interval.type22[pos]
  } else if (interval.type %in% interval.type23) {
    interval.type <- interval.type23[pos]
  } else if (interval.type %in% interval.type24) {
    interval.type <- interval.type24[pos]
  } else if (interval.type %in% interval.type25) {
    interval.type <- interval.type25[pos]
  } else if (interval.type %in% interval.type31) {
    interval.type <- interval.type31[pos]
  } else if (interval.type %in% interval.type32) {
    interval.type <- interval.type32[pos]
  } else if (interval.type %in% interval.type33) {
    interval.type <- interval.type33[pos]
  } else if (interval.type %in% interval.type34) {
    interval.type <- interval.type34[pos]
  } else if (interval.type %in% interval.type35) {
    interval.type <- interval.type35[pos]
  } else if (interval.type %in% interval.type36) {
    interval.type <- interval.type36[pos]
  } else if (interval.type %in% interval.type37) {
    interval.type <- interval.type37[pos]
  } else if (interval.type %in% interval.type41) {
    interval.type <- interval.type41[pos]
  } else if (interval.type %in% interval.type42) {
    interval.type <- interval.type42[pos]
  } else if (interval.type %in% interval.type43) {
    interval.type <- interval.type43[pos]
  } else if (interval.type %in% interval.type44) {
    interval.type <- interval.type44[pos]
  } else {
    stop("Please provide a correct type of confidence interval. Please see the help page. However, if we were to recommend one, we recommend the 'boot' method.")
  }
  interval.type
}

