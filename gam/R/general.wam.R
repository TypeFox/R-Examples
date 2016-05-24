"general.wam" <-
  function(x, y, w, s, which, smooth.frame, maxit = 30, tol = 1e-7, trace = FALSE,
           se = TRUE, ...)
{
  if(inherits(smooth.frame, "data.frame")) {
    data <- smooth.frame
### Note; the lev component of the smooths is the diagonal hat matrix elements
### for the NONLINEAR part of the fit.
###The smoother can return both the linear and nonlinear parts, although only
### the nonlinear part is strictly necessary. 
###
    oldClass(data) <- NULL
    names.calls <- names(which)
    smooth.calls <- lapply(data[names.calls], attr, "call")
    names(smooth.calls) <- names.calls
    smooth.frame <- list(data = data, smooth.calls = smooth.calls)
  }
  else {
    data <- smooth.frame$data
    smooth.calls <- smooth.frame$smooth.calls
  }
  names.calls <- names(smooth.calls)
  y <- as.vector(y)
  residuals <- as.vector(y - s %*% rep(1., ncol(s)))
  n <- length(y)
  fit <- list(fitted.values = 0.)
  rss <- weighted.mean(residuals^2., w)
  rssold <- rss * 10.
  nit <- 0.
  df <- rep(NA, length(which))
  var <- s
  if(trace)
    cat("\nWAM   iter   rss/n     term\n")
  ndig <-  - log10(tol) + 1.
  RATIO <- tol + 1.
  while(RATIO > tol & nit < maxit) {
    rssold <- rss
    nit <- nit + 1.
    z <- residuals + fit$fitted.values
    fit <- lm.wfit(x, z, w, method = "qr", singular.ok = TRUE,
                   ...)
    residuals <- fit$residuals
    rss <- weighted.mean(residuals^2., w)
    if(trace)
      cat("\n         ", nit, "   ", format(round(rss, ndig)),
          "  Parametric -- lm.wfit\n", sep = "")
    deltaf <- 0.
    for(j in seq(names.calls)) {
      old <- s[, j]
      z <- residuals + s[, j]
      fit.call <- eval(smooth.calls[[j]])
      residuals <- as.double(fit.call$residuals)
      if(length(residuals) != n)
        stop(paste(names.calls[j], 
                   "returns a vector of the wrong length")
             )
      s[, j] <- z - residuals
      deltaf <- deltaf + weighted.mean((s[, j] - old)^2.,
                                       w)
      rss <- weighted.mean(residuals^2., w)
      if(trace) {
        cat("         ", nit, "   ", format(round(
                                                  rss, ndig)), "  Nonparametric -- ",
            names.calls[j], "\n", sep = "")
      }
      df[j] <- fit.call$nl.df
      if(se)
        var[, j] <- fit.call$var
    }
    RATIO <- sqrt(deltaf/sum(w * apply(s, 1., sum)^2.))
    if(trace)
      cat("Relative change in functions:", format(round(
                                                        RATIO, ndig)), "\n")
  }
  if((nit == maxit) & maxit > 1.)
    warning(paste("general.wam convergence not obtained in ", maxit,
                  " iterations"))
  names(df) <- names.calls
  if(trace)
    cat("\n")
  fit$fitted.values <- y - residuals
  rl <- c(fit, list(smooth = s, nl.df = df))
  rl$df.residual <- rl$df.residual - sum(df)
  if(se)
    rl <- c(rl, list(var = var))
  c(list(smooth.frame = smooth.frame), rl)
}
