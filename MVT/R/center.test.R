center.test <-
function(object, center, test = "LRT")
{
  ## local functions
  restrictedCenter <-
  function(z, dims, settings, center, Scatter, control)
  {
    ctrl <- unlist(control)
    ctrl <- c(ctrl, 0)
    n <- dims[1]
    p <- dims[2]
    o <- .C("restricted_center",
            z = as.double(t(z)),
            dims = as.integer(dims),
            settings = as.double(settings),
            center = as.double(center),
            Scatter = as.double(Scatter),
            distances = double(n),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            control = as.double(ctrl))
    o$Scatter <- matrix(o$Scatter, ncol = p)
    o$eta <- o$settings[2]
    o$numIter <- o$control[4]
    o
  }
  restrictedScore <-
  function(z, dims, weights, center)
  { # compute the unscaled score
    n <- dims[1]
    p <- dims[2]
    y <- z - matrix(rep(center, times = n), ncol = p, byrow = TRUE)
    y <- weights * y
    s <- apply(y, 2, sum)
    s
  }
  
  ## extract object info
  fit <- object
  if (!inherits(fit, "studentFit"))
    stop("Use only with 'studentFit' objects")
  n <- fit$dims[1]
  p <- fit$dims[2]
  if (length(center) != p)
    stop("center must have", p, "elements")
  names(center) <- names(fit$center)
  
  ## restricted parameter estimation
  f0 <- restrictedCenter(fit$x, fit$dims, fit$settings, center, fit$Scatter, fit$control)
  null.fit <- list(center = f0$center, Scatter = f0$Scatter, eta = f0$eta, distances = f0$distances,
                   weights = f0$weights, logLik = f0$logLik, numIter = f0$numIter)

  switch(test,
        hotelling = {
          dif <- colMeans(fit$x) - center
          R <- chol(cov(fit$x))
          dif <- solve(t(R), dif)
          stat <- n * sum(dif^2)
          names(stat) <- "Hotelling"
          method <- "Hotelling's T-squared test"
        },
        LRT = {
          stat <- 2. * (fit$logLik - f0$logLik)
          names(stat) <- "LRT"
          method <- "Likelihood ratio test"
        },
        Wald = {
          dif <- fit$center - center
          aCov <- fisher.info(fit)[1:p,1:p]
          R <- chol(aCov)
          dif <- R %*% dif
          stat <- n * sum(dif^2)
          names(stat) <- "Wald"
          method <- "Wald test"
        },
        score = {
          s <- restrictedScore(fit$x, fit$dims, f0$weights, center)
          s <- solve(f0$Scatter, s)
          aCov <- fisher.info(f0)[1:p,1:p]
          R <- chol(aCov)
          s <- solve(t(R), s)
          stat <- sum(s^2) / n
          names(stat)<-"Score"
          method <- "Score test"
        },
        gradient = {
          dif <- fit$center - center
          dif <- solve(f0$Scatter, dif)
          s <- restrictedScore(fit$x, fit$dims, f0$weights, center)
          stat <- sum(s * dif)
          names(stat) <- "Gradient"
          method <- "Gradient test"
        },
        stop(paste("unimplemented test:", test))
        )
  pval <- 1 - pchisq(stat, df = p)
  
  ## output object
  z <- list(statistic = stat, parameter = p, p.value = pval, estimate = fit$center,
            null.value = center, method = method, null.fit = null.fit)
  if (is.null(fit$call$data))
    z$data <- fit$call$x
  else
    z$data <- fit$call$data
  z$family <- fit$family
  class(z) <- "center.test"
  z
}

print.center.test <- function(x, digits = 4, ...)
{
  cat("\n")
  cat(paste(x$method, "for the center parameter", sep = " "), "\n")
  cat("\n")
  cat("data:", x$data, "\n")
  out <- character()
  out <- c(out, paste(names(x$statistic), "statistic =",
                      format(round(x$statistic, digits = digits))))
  out <- c(out, paste("df =", x$parameter))
  out <- c(out, paste("p-value =", format(round(x$p.value, digits = digits))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: true center parameter is not equal to:\n")
  print(x$null.value, ...)
  cat("sample estimates:\n")
  print(x$estimate, ...)
  invisible(x)
}
