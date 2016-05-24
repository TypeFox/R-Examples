print.spatgev <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("      Model: Spatial GEV model\n")
  cat("   Deviance:", x$deviance, "\n")
  cat("        TIC:", TIC(x), "\n\n")

  param <- x$fitted.values
  loc.idx <- which(substr(names(param), 1, 3) == "loc")
  scale.idx <- which(substr(names(param), 1, 5) == "scale")
  shape.idx <- which(substr(names(param), 1, 5) == "shape")

  cat("    Location Parameters:\n")
  print.default(format(param[loc.idx], digits = digits), print.gap = 2,
                quote = FALSE)
  cat("       Scale Parameters:\n")
  print.default(format(param[scale.idx], digits = digits), print.gap = 2,
                quote = FALSE)
  cat("       Shape Parameters:\n")
  print.default(format(param[shape.idx], digits = digits), print.gap = 2,
                quote = FALSE)

  loc.idx <- which(substr(names(param), 1, 12) == "tempCoeffLoc")
  scale.idx <- which(substr(names(param), 1, 14) == "tempCoeffScale")
  shape.idx <- which(substr(names(param), 1, 14) == "tempCoeffShape")

  if ((length(loc.idx) + length(scale.idx) + length(shape.idx)) > 0){
    cat("\nTemporal Parameters\n")

    if (length(loc.idx) > 0){
      cat("    Location Parameters:\n")
      print.default(format(param[loc.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }

    if (length(scale.idx)> 0){
      cat("       Scale Parameters:\n")
      print.default(format(param[scale.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }

    if (length(shape.idx)>0){
      cat("       Shape Parameters:\n")

      print.default(format(param[shape.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }
  }

  if(!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

print.maxstab <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("        Estimator:", x$est, "\n")
  cat("            Model:", x$model, "\n")
  cat("         Weighted:", x$weighted, "\n")
  if (x$est == 'MPLE'){
    cat("   Pair. Deviance:", x$deviance, "\n")
    cat("              TIC:", TIC(x), "\n")
  }
  if (x$est == "Least Squares")
    cat("  Objective Value:", x$opt.value, "\n")

  if ((x$model == "Schlather") || (x$model == "Geometric") || (x$model == "Brown-Resnick") ||
      (x$model == "Extremal-t")){

    if (x$cov.mod == "emp")
      cov.mod <- "Empirical"

    if (x$cov.mod == "whitmat")
      cov.mod <- "Whittle-Matern"

    if (x$cov.mod == "powexp")
      cov.mod <- "Powered Exponential"

    if (x$cov.mod == "cauchy")
      cov.mod <- "Cauchy"

    if (x$cov.mod == "caugen")
      cov.mod <- "Generalized Cauchy"

    if (x$cov.mod == "bessel")
      cov.mod <- "Bessel"

    if (x$cov.mod == "brown")
      cov.mod <- "Fractional Brownian Motion"

    cat("Covariance Family:", cov.mod, "\n")

    cat("\nEstimates\n")
    cat("  Marginal Parameters:\n")

    if (x$fit.marge){
      idx <- which(names(x$fitted.values) == "alpha")
      idx <- c(idx, which(names(x$fitted.values) == "sigma2"))
      idx <- c(idx, which(names(x$fitted.values) == "nugget"))
      idx <- c(idx, which(names(x$fitted.values) == "range"))
      idx <- c(idx, which(names(x$fitted.values) == "smooth"))
      idx <- c(idx, which(names(x$fitted.values) == "smooth2"))
      idx <- c(idx, which(names(x$fitted.values) == "DoF"))

      margin.param <- x$fitted.values[-idx]
      loc.idx <- which(substr(names(margin.param), 1, 3) == "loc")
      scale.idx <- which(substr(names(margin.param), 1, 5) == "scale")
      shape.idx <- which(substr(names(margin.param), 1, 5) == "shape")

      cat("    Location Parameters:\n")
      print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
      cat("       Scale Parameters:\n")
      print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
      cat("       Shape Parameters:\n")
      print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2,
                    quote = FALSE)

      loc.idx <- which(substr(names(margin.param), 1, 12) == "tempCoeffLoc")
      scale.idx <- which(substr(names(margin.param), 1, 14) == "tempCoeffScale")
      shape.idx <- which(substr(names(margin.param), 1, 14) == "tempCoeffShape")

      if (length(loc.idx) > 0){
        cat("Temporal Location Parameters:\n")
        print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2,
                      quote = FALSE)
      }

      if (length(scale.idx)> 0){
        cat("Temporal Scale Parameters:\n")
        print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2,
                      quote = FALSE)
      }

      if (length(shape.idx)>0){
        cat("Temporal Shape Parameters:\n")

        print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2,
                      quote = FALSE)
      }

      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values[idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }

    else{
      cat("  Assuming unit Frechet.\n\n")
      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values, digits = digits), print.gap = 2,
                    quote = FALSE)
    }
  }

  else{
    cat("Covariance Family:", x$cov.mod, "\n")

    cat("\nEstimates\n")
    cat("  Marginal Parameters:\n")

    if (x$fit.marge){
      idx <- which(substr(names(x$fitted.values), 1, 3) == "cov")

      margin.param <- x$fitted.values[-idx]
      loc.idx <- which(substr(names(margin.param), 1, 3) == "loc")
      scale.idx <- which(substr(names(margin.param), 1, 5) == "scale")
      shape.idx <- which(substr(names(margin.param), 1, 5) == "shape")

      cat("    Location Parameters:\n")
      print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
      cat("       Scale Parameters:\n")
      print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
      cat("       Shape Parameters:\n")
      print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2,
                    quote = FALSE)

      loc.idx <- which(substr(names(margin.param), 1, 12) == "tempCoeffLoc")
      scale.idx <- which(substr(names(margin.param), 1, 14) == "tempCoeffScale")
      shape.idx <- which(substr(names(margin.param), 1, 14) == "tempCoeffShape")

      if (length(loc.idx) > 0){
        cat("Temporal Location Parameters:\n")
        print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2,
                      quote = FALSE)
      }

      if (length(scale.idx)> 0){
        cat("Temporal Scale Parameters:\n")
        print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2,
                      quote = FALSE)
      }

      if (length(shape.idx)>0){
        cat("Temporal Shape Parameters:\n")

        print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2,
                      quote = FALSE)
      }

      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values[idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }

    else{
      cat("  Not estimated.\n")
      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values, digits = digits), print.gap = 2,
                    quote = FALSE)
    }
  }

  if(!is.null(x$std.err)) {
      cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

logLik.maxstab <- function(object, ...){
  llk <- object$logLik
  attr(llk, "df") <- length(fitted(object))
  class(llk) <- "logLik"
  return(llk)
}

logLik.copula <- function(object, ...){
  llk <- object$logLik
  attr(llk, "df") <- length(fitted(object))
  class(llk) <- "logLik"
  return(llk)
}

profile2d <- function(fitted, ...){
  UseMethod("profile2d")
}

print.pspline <- function(x, ...){
  cat("Call:\n")
  print(x$call)

  cat("\n  Rank:", x$rank, "\t(G)CV Score:", round(x$cv, 3),
      "\n")
  cat("Degree:", x$degree, "\t Penalty: ",
      round(x$penalty, 3), "\n")
  cat("\n     Degree of freedom:", round(x$df, 3), "\n")
  cat("Res. Degree of freedom:", round(x$res.df, 3), "\n")
}

TIC <- function(object, ..., k = 2){
  UseMethod("TIC")
}

print.latent <- function(x, digits = max(3, getOption("digits") - 3), ...,
                         level = 0.95){

  if ((level > 1) || (level < 0))
    stop("'level' must lie in [0, 1]")

  alpha <- 0.5 * (1 - level)

  cat("Effective length:", nrow(x$chain.loc), "\n")
  cat("         Burn-in:", x$burn.in, "\n")
  cat("        Thinning:", x$thin, "\n")
  cat("   Effective NoP:", x$eNoP, "\n")
  cat("             DIC:", x$DIC, "\n\n")


  cat("  Regression Parameters:\n")
  cat("      Location Parameters:\n")
  chain <- x$chain.loc
  loc.idx <- which(substr(colnames(chain), 1, 2) == "lm")
  post.mean <- colMeans(chain[,loc.idx, drop=FALSE])
  ci.lower <- apply(chain[,loc.idx, drop=FALSE], 2, quantile, alpha)
  ci.upper <- apply(chain[,loc.idx, drop=FALSE], 2, quantile, 1-alpha)
  dummy <- rbind(ci.lower = ci.lower, post.mean = post.mean,
                 ci.upper = ci.upper)
  print.default(format(dummy, digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")
  cat("         Scale Parameters:\n")
  chain <- x$chain.scale
  scale.idx <- which(substr(colnames(chain), 1, 2) == "lm")
  post.mean <- colMeans(chain[,scale.idx, drop=FALSE])
  ci.lower <- apply(chain[,scale.idx, drop=FALSE], 2, quantile, alpha)
  ci.upper <- apply(chain[,scale.idx, drop=FALSE], 2, quantile, 1-alpha)
  dummy <- rbind(ci.lower = ci.lower, post.mean = post.mean,
                 ci.upper = ci.upper)
  print.default(format(dummy, digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")
  cat("         Shape Parameters:\n")
  chain <- x$chain.shape
  shape.idx <- which(substr(colnames(chain), 1, 2) == "lm")
  post.mean <- colMeans(chain[,shape.idx, drop=FALSE])
  ci.lower <- apply(chain[,shape.idx, drop=FALSE], 2, quantile, alpha)
  ci.upper <- apply(chain[,shape.idx, drop=FALSE], 2, quantile, 1-alpha)
  dummy <- rbind(ci.lower = ci.lower, post.mean = post.mean,
                 ci.upper = ci.upper)
  print.default(format(dummy, digits = digits), print.gap = 2, quote = FALSE)

  cat("\n\n")
  cat("  Latent Parameters:\n")
  cat("      Location Parameters:\n")
  cat("        Covariance family:", x$cov.mod[1],"\n")
  chain <- x$chain.loc
  loc.idx <- which(colnames(chain) %in% c("sill", "range", "smooth"))
  post.mean <- colMeans(chain[,loc.idx, drop=FALSE])
  ci.lower <- apply(chain[,loc.idx, drop=FALSE], 2, quantile, alpha)
  ci.upper <- apply(chain[,loc.idx, drop=FALSE], 2, quantile, 1-alpha)
  dummy <- rbind(ci.lower = ci.lower, post.mean = post.mean,
                 ci.upper = ci.upper)
  print.default(format(dummy, digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")
  cat("      Scale Parameters:\n")
  cat("     Covariance family:", x$cov.mod[2],"\n")
  chain <- x$chain.scale
  scale.idx <- which(colnames(chain) %in% c("sill", "range", "smooth"))
  post.mean <- colMeans(chain[,scale.idx, drop=FALSE])
  ci.lower <- apply(chain[,scale.idx, drop=FALSE], 2, quantile, alpha)
  ci.upper <- apply(chain[,scale.idx, drop=FALSE], 2, quantile, 1-alpha)
  dummy <- rbind(ci.lower = ci.lower, post.mean = post.mean,
                 ci.upper = ci.upper)
  print.default(format(dummy, digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")
  cat("      Shape Parameters:\n")
  cat("     Covariance family:", x$cov.mod[3],"\n")
  chain <- x$chain.shape
  shape.idx <- which(colnames(chain) %in% c("sill", "range", "smooth"))
  post.mean <- colMeans(chain[,shape.idx, drop=FALSE])
  ci.lower <- apply(chain[,shape.idx, drop=FALSE], 2, quantile, alpha)
  ci.upper <- apply(chain[,shape.idx, drop=FALSE], 2, quantile, 1-alpha)
  dummy <- rbind(ci.lower = ci.lower, post.mean = post.mean,
                 ci.upper = ci.upper)
  print.default(format(dummy, digits = digits), print.gap = 2, quote = FALSE)
}

print.copula <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("           Copula:", x$copula, "\n")
  cat("         Deviance:", x$deviance, "\n")
  cat("              AIC:", AIC(x), "\n")

  if (x$cov.mod == "whitmat")
    cov.mod <- "Whittle-Matern"

  if (x$cov.mod == "powexp")
    cov.mod <- "Powered Exponential"

  if (x$cov.mod == "cauchy")
    cov.mod <- "Cauchy"

  if (x$cov.mod == "caugen")
    cov.mod <- "Generalized Cauchy"

  if (x$cov.mod == "bessel")
    cov.mod <- "Bessel"

  cat("Covariance Family:", cov.mod, "\n")

  idx <- which(names(x$fitted.values) == "DoF")
  idx <- c(idx, which(names(x$fitted.values) == "nugget"))
  idx <- c(idx, which(names(x$fitted.values) == "range"))
  idx <- c(idx, which(names(x$fitted.values) == "smooth"))
  idx <- c(idx, which(names(x$fitted.values) == "smooth2"))


  cat("\nEstimates\n")
  cat("  Marginal Parameters:\n")

  if (x$fit.marge){
    margin.param <- x$fitted.values[-idx]
    loc.idx <- which(substr(names(margin.param), 1, 3) == "loc")
    scale.idx <- which(substr(names(margin.param), 1, 5) == "scale")
    shape.idx <- which(substr(names(margin.param), 1, 5) == "shape")

    cat("    Location Parameters:\n")
    print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("       Scale Parameters:\n")
    print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("       Shape Parameters:\n")
    print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2,
                  quote = FALSE)

    loc.idx <- which(substr(names(margin.param), 1, 12) == "tempCoeffLoc")
    scale.idx <- which(substr(names(margin.param), 1, 14) == "tempCoeffScale")
    shape.idx <- which(substr(names(margin.param), 1, 14) == "tempCoeffShape")

    if (length(loc.idx) > 0){
      cat("Temporal Location Parameters:\n")
      print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }

    if (length(scale.idx)> 0){
      cat("Temporal Scale Parameters:\n")
      print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }

    if (length(shape.idx)>0){
      cat("Temporal Shape Parameters:\n")

      print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2,
                    quote = FALSE)
    }
  }

  else
    cat("  Assuming unit Frechet.\n\n")

  cat("  Dependence Parameters:\n")
  print.default(format(x$fitted.values[idx], digits = digits), print.gap = 2,
                quote = FALSE)

  if(!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

