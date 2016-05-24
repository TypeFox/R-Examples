
coef.stsmFit <- function(object, ...) {
  rescale <- !is.null(object$model@cpar)
  pars <- c(get.pars(object$model, rescale), get.cpar(object$model, rescale))
  pars[sort(names(pars))]
}

print.stsmFit <- function(x, digits = max(3L, getOption("digits") - 3L), 
  vcov.type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), ...)
{  
  cat("Call:", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
  cat("Parameter estimates:\n")
  rescale <- !is.null(x$model@cpar)
  pars <- c(get.pars(x$model, rescale), get.cpar(x$model, rescale))
  ord <- sort(names(pars))

  # if default method for each method is kept

  if (length(vcov.type) > 1)
  {
    res <- rbind(pars[ord], x$std.errors[ord])
    rownames(res) <- c("Estimate", "Std. error")
    print.default(res, print.gap = 2L, digits = digits, ...)

    cat("\nLog-likelihood:", 
      round(x$loglik, digits = digits), "\n")
    #print.default(x$loglik, print.gap = 2L, digits = digits, ...)

    if (grepl("optim", x$call[1])) {
      cat("Convergence code:", as.numeric(x$convergence), "\n")
      if (!is.null(x$message) && x$message != "")
        cat(x$message, "\n")
    } else
      cat("Convergence:", x$convergence, "\n")

    if (grepl("optim", x$call[1])) {
      cat("Number of function calls:", x$iter[1], "\n")
    } else
      cat("Number of iterations:", x$iter, "\n")

    cat("Variance-covariance matrix:", x$vcov.type, "\n")

    #invisible(x)
    return()
  }

  # if "vcov.type" is specified by the user

  vcov.type <- match.arg(vcov.type)

  #NOTE do not pass 'vcov.type[1]' since method 'vcov()' uses length(vcov.type)
  #to see whether a choice was made or not and then set default value
  if (x$model@model != "trend+ar2")
  {
##FIXME vcov is computed whenever print is called
    res <- rbind(pars[ord], sqrt(diag(vcov(x, type = vcov.type)))[ord])
    rownames(res) <- c("Estimate", "Std. error")
    print.default(res, print.gap = 2L, digits = digits, ...)
  } else {
    res <- rbind(pars[ord])
    rownames(res) <- c("Estimate")
    print.default(res, print.gap = 2L, digits = digits, ...)
  }

##FIXME TODO std.error
  if (!is.null(x$xreg))
  {
    cat("\nCoefficients of regressor variables:\n")
    if (is.null(x$xreg$stde))
    {
      if (vcov.type == "optimHessian") {
        # this requires having run maxlik.td.optim() with option hessian=TRUE
        stde <- sqrt(diag(vcov(x, type = "optimHessian"))[names(x$xreg$coef)])
        res <- rbind(x$xreg$coef, stde)
        rownames(res) <- c("Estimate", "Std. error")
      } else {
        res <- rbind(x$xreg$coef)
        rownames(res) <- c("Estimate")
      }
    } else {
      res <- rbind(x$xreg$coef, x$xreg$stde)
      rownames(res) <- c("Estimate", "Std. error")
    }
    print.default(res, print.gap = 2L, digits = digits, ...)
  }

  cat("\nLog-likelihood:", 
    round(x$loglik, digits = digits), "\n")
  #print.default(x$loglik, print.gap = 2L, digits = digits, ...)

  cat("Convergence:", as.numeric(x$convergence), "\n")
  if (!is.null(x$message) && x$message != "")
    cat(x$message, "\n")

  if (grepl("optim", x$call[1])) {
    cat("Number of function calls:", x$iter[1], "\n")
  } else
    cat("Number of iterations:", x$iter, "\n")

##FIXME print method used for std.errors

  invisible(x)
}

fitted.stsm <- function(object, std.rediduals = TRUE, 
  version = c("KFKSDS", "stats"), ...)
{
  version <- match.arg(version)[1]

  y <- object@y
  id <- which(diag(object@ss$Q) != "0")
  ss <- char2numeric(object, ...)

  switch(version,

    "KFKSDS" = {
      kf <- KFKSDS::KF(y, ss)

##FIXME return matrix in "llm" model for kf$a.upd in KFKSDS, arrange also kf$P.upd
      if (is.null(dim(kf$a.upd))) {
        states <- cbind(kf$a.upd[id])
      } else
        states <- kf$a.upd[,id]
      P <- matrix(nrow = length(y), ncol = length(id))
      if (is.null(dim(kf$P.upd))) {
        for (i in seq_along(y))
          P[i,] <- sqrt(kf$P.upd[i])
      } else
        for (i in seq_along(y))
          P[i,] <- sqrt(diag(kf$P.upd[id,id,i]))
      P <- ts(P, start = start(y)[1L], frequency = frequency(y))

      if (std.rediduals) {
        resid <- kf$v / sqrt(kf$f)
      } else
        resid <- kf$v
    },

    "stats" = { # based on stats::StructTS
      mod <- list(Z = ss$Z, a = ss$a0, P = ss$P0, T =  ss$T, 
        V =  ss$Q, h =  ss$H, Pn = ss$P0)
      z <- stats::KalmanRun(y, mod, nit = -1, fast = TRUE)

      resid <- z$resid
      states <- z$states
      P <- NULL
    }
  )

  resid <- ts(resid, start = start(y)[1L], frequency = frequency(y))
  states <- ts(states, start = start(y)[1L], frequency = frequency(y))

  nms <- grep("^var\\d{1,2}$", names(object@pars), value = TRUE)
  if ("var1" %in% nms)
    nms <- nms[-which(nms == "var1")]
  colnames(states) <- nms

  # return here the states and also the residuals 
  # to avoid running the filter two separate times;
  # computing the states involves computing the residuals and vice versa,
  # in 'StructTS' 'fitted' and 'residuals' take the objects in class 'StructTS'
  # here those computations are a separate task and 
  # they are not done within 'maxlik.fd.scoring'

  res <- list(states = states, states.se = P, residuals = resid)
  class(res) <- "stsmComponents"
  res
}

fitted.stsmFit <- function(object, std.rediduals = TRUE, 
  version = c("KFKSDS", "stats"), ...)
{
  fitted(object$model, std.rediduals, version, ...)
}

residuals.stsmFit <- function(object, standardised = FALSE, 
  version = c("KFKSDS", "stats"), ...)
{
  # the computations are the same as those for obtaining the states,
  # computing the states involves computing the residuals and vice versa,
  # so method 'fitted' is reused;
  # in 'StructTS' the residuals are taking the output from the objects 
  # available in class 'StructTS';
  # here those computations are a separate task and 
  # they are not done within 'maxlik.fd.scoring'

  fitted(object, standardised, version, ...)$residuals
}

plot.stsmComponents <- function(x, ...)
{
  states <- x$states
  se <- x$states.se

  if (is.ts(states)) {
    xpol <- c(time(states), rev(time(states)))
  } else
    xpol <- c(seq_len(nrow(states)), seq.int(nrow(states), 1))

  sse1 <- states + 1.96 * se
  sse2 <- states - 1.96 * se

  oldpar <- par(mfrow = c(ncol(states), 1))
  on.exit(par(oldpar))
  for (i in seq_len(ncol(states)))
  {
    plot(ts.union(states[,i], sse1[,i], sse2[,i]),
      ylab = "", plot.type = "single", type = "n", ...)
    #lines(sse1[,i], col = 2, lty = 2)
    #lines(sse2[,i], col = 2, lty = 2)
    polygon(x = xpol, y = c(sse2[,i], rev(sse1[,i])), 
      col = "lightgray", border = NA)
    lines(states[,i])
  }
  #invisible(x)
}

##FIXME see function in package KFKSDS

predict.stsm <- function(object, n.ahead = 1L, se.fit = TRUE, 
  version = c("KFKSDS", "stats"), ...)
{
  version <- match.arg(version)[1]

  y <- object@y
  nobs <- length(y)
  xtsp <- tsp(y)

  ss <- char2numeric(object, ...)

  switch(version,

    "KFKSDS" = {

      kf <- KFKSDS::KF(y, ss)
      a <- ss$T %*% kf$a.upd[nobs,]
      P <- ss$T %*% kf$P.upd[,,nobs] %*% t(ss$T) + ss$Q

      if (se.fit)
      {
        ypred <- se <- rep(NA, n.ahead)
        for (i in seq_len(n.ahead))
        {
          ypred[i] <- sum(ss$Z * drop(a))
          a <- ss$T%*% a

          se[i] <- sqrt(drop(ss$H + ss$Z %*% P %*% t(ss$Z)))
          P <- ss$T %*% P %*% t(ss$T) + ss$Q
        }

        if (is.ts(y))
        {
          ypred <- ts(ypred, start = xtsp[2L] + 1/xtsp[3L], frequency = xtsp[3L])
          se <- ts(se, start = xtsp[2L] + 1/xtsp[3L], frequency = xtsp[3L])
        }

        res <- list(pred = ypred, se = se, y = y)

      } else 
      {
        ypred <- rep(NA, n.ahead)
        for (i in seq_len(n.ahead))
        {
          ypred[i] <- sum(ss$Z * drop(a))
          a <- ss$T%*% a
        }

        if (is.ts(y))
        {
          ypred <- ts(ypred, start = xtsp[2L] + 1/xtsp[3L], frequency = xtsp[3L])
        }

        res <- list(pred = ypred, se = NULL, y = y)
      }
    },

    "stats" = { # based on stats::KalmanForecast

      mod <- list(Z = ss$Z, a = ss$a0, P = ss$P0, T =  ss$T, 
      V =  ss$Q, h =  ss$H, Pn = ss$P0)
      z <- stats::KalmanRun(y, mod, nit = 0, fast = TRUE)
      # 'mod' is modified by 'KalmanRun', in particular mod$a 
      # contains the last state prediction and mod$P its covariance matrix
      z <- stats::KalmanForecast(n.ahead, mod, fast = TRUE)
      pred <- ts(z[[1L]], start = xtsp[2L] + 1/xtsp[3L], frequency = xtsp[3L])

      if (se.fit)
      {
        se <- ts(sqrt(z[[2L]]), start = xtsp[2L] + 1/xtsp[3L], frequency = xtsp[3L])
        res <- list(pred = pred, se = se, y = y)
      }
      else 
        res <- list(pred = pred, se = NULL, y = y)
    }
  )

  class(res) <- "stsmPredict"
  res
}

predict.stsmFit <- function(object, n.ahead = 1L, se.fit = TRUE, 
  version = c("KFKSDS", "stats"), ...)
{
  predict(object$model, n.ahead = n.ahead, se.fit = se.fit, version = version, ...)
}

plot.stsmPredict <- function(x, ...)
{
  pred <- x$pred 
  se <- x$se
  y <- x$y

  predse1 <- pred + 1.96 * se
  predse2 <- pred - 1.96 * se

  plot(ts.union(y, predse1, predse2), 
    ylab = "", plot.type = "single", type = "n", ...)
  lines(y, col = 1)
  lines(pred, col = 4)
  lines(predse1, col = 2, lty = 2)
  lines(predse2, col = 2, lty = 2)
  #invisible(x)
}

tsSmooth.stsm <- function(object, version = c("KFKSDS", "stats"), ...)
{
  version <- match.arg(version)[1]

  y <- object@y
  id <- which(diag(object@ss$Q) != "0")
  ss <- char2numeric(object, ...)

  switch(version,

    "KFKSDS" = {
      
      kf <- KFKSDS::KF(y, ss)
      ks <- KFKSDS::KS(y, ss, kf)
      #ks <- KSDS(y, ss, kf)
      #ks <- KalmanSmoother()

      states <- ks$ahat[,id]
      sse <- matrix(nrow = nrow(states), ncol = length(id))
      for (i in seq_len(nrow(states)))
      {
        sse[i,] <- sqrt(diag(ks$varahat[id,id,i]))
      }

    },

    "stats" = { # based on stats::StructTS

      mod <- list(Z = ss$Z, a = ss$a0, P = ss$P0, T =  ss$T, 
        V =  ss$Q, h =  ss$H, Pn = ss$P0)
      res <- KalmanSmooth(object$model@y, mod, nit = -1)

      states <- res$smooth[,id]
      states <- ts(states, start = start(y)[1L], frequency = frequency(y))

      sse <- matrix(nrow = nrow(states), ncol = length(id))
      for (i in seq_len(nrow(states)))
      {
        sse[i,] <- sqrt(diag(res$var[i,id,id]))
      }      
    }
  )

  sse <- ts(sse, start = start(y)[1L], frequency = frequency(y))

  nms <- grep("^var\\d{1,2}$", names(object@pars), value = TRUE)
  if ("var1" %in% nms)
    nms <- nms[-which(nms == "var1")]
  colnames(states) <- colnames(sse) <- nms

  res <- list(states = states, sse = sse)
  class(res) <- "stsmSmooth"
  res
}

tsSmooth.stsmFit <- function(object, version = c("KFKSDS", "stats"), ...)
{
  tsSmooth(object$model, version = version, ...)
}

plot.stsmSmooth <- function(x, ...)
{
  states <- x$states
  se <- x$sse

  if (is.ts(states)) {
    xpol <- c(time(states), rev(time(states)))
  } else
    xpol <- c(seq_len(nrow(states)), seq.int(nrow(states), 1))

  sse1 <- states + 1.96 * se
  sse2 <- states - 1.96 * se

  oldpar <- par(mfrow = c(ncol(states), 1), ...)
  on.exit(par(oldpar))
  for (i in seq_len(ncol(states)))
  {
    plot(ts.union(states[,i], sse1[,i], sse2[,i]),
      ylab = "", plot.type = "single", type = "n", ...)
    #lines(sse1[,i], col = 2, lty = 2)
    #lines(sse2[,i], col = 2, lty = 2)
    polygon(x = xpol, y = c(sse2[,i], rev(sse1[,i])), 
      col = "lightgray", border = NA)
    lines(states[,i])
  }
  #invisible(x)
}

tsdiag.stsm <- function(object, gof.lag = 10L, ...)
{
#adapted from 'tsdiag.StructTS'

    # plot standardized residuals, acf of residuals, Ljung-Box p-values
    oldpar <- par(mfrow = c(3, 1))
    on.exit(par(oldpar))
    stdres <- residuals(object)
    plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
    abline(h = 0.0)    
    acf(stdres, plot = TRUE, main = "ACF of Residuals", na.action = na.pass)
    nlag <- gof.lag
    pval <- numeric(nlag)
    for(i in 1L:nlag)
      pval[i] <- Box.test(stdres, i, type = "Ljung-Box")$p.value
    plot(1L:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0,1),
         main = "p values for Ljung-Box statistic")
    abline(h = 0.05, lty = 2L, col = "blue")
}

tsdiag.stsmFit <- function(object, gof.lag = 10L, ...)
{
#adapted from 'tsdiag.StructTS'

    # plot standardized residuals, acf of residuals, Ljung-Box p-values
    oldpar <- par(mfrow = c(3, 1))
    on.exit(par(oldpar))
    stdres <- residuals(object)
    plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
    abline(h = 0.0)    
    acf(stdres, plot = TRUE, main = "ACF of Residuals", na.action = na.pass)
    nlag <- gof.lag
    pval <- numeric(nlag)
    for(i in 1L:nlag)
      pval[i] <- Box.test(stdres, i, type = "Ljung-Box")$p.value
    plot(1L:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0,1),
         main = "p values for Ljung-Box statistic")
    abline(h = 0.05, lty = 2L, col = "blue")
}
