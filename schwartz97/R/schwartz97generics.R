### <======================================================================>
resid.schwartz2f.fit <- function(object, data, ttm, type = c("filter", "filter.std", "real"))
{
  type <- match.arg(type)
  if(type == "filter"){
    vt <- t(filter.schwartz2f(data, ttm, object)$fkf.obj$vt)
    dimnames(vt) <- dimnames(data)
    return(vt)
  }else if(type == "filter.std"){
    filter.obj <- filter.schwartz2f(data, ttm, object)$fkf.obj
    resid.std <- t(sapply(1:ncol(filter.obj$vt), function(i, Ft, vt){if(any(is.na(vt[,i])))
                                                                       return(matrix(NA, ncol = 1, nrow = nrow(vt)))
                                                                     else
                                                                       return(solve(t(chol(Ft[,,i]))) %*% vt[,i])
                                                                     },
                          Ft = filter.obj$Ft, vt = filter.obj$vt))
    dimnames(resid.std) <- dimnames(data)
    return(resid.std)
  }else if(type == "real"){
    resid.real <- data - fitted(object, data, ttm)
    return(resid.real)
  }
}
### <---------------------------------------------------------------------->
setMethod("resid", signature(object = "schwartz2f.fit"), resid.schwartz2f.fit)
### <---------------------------------------------------------------------->

### <======================================================================>
fitted.schwartz2f.fit <- function(object, data, ttm)
{

  state <- filter.schwartz2f(data, ttm, object)
  coefs <- coef(object)

  pricefutures.wrapper <- function(x, mu, sigmaS, kappa, alpha, sigmaE, rho, r, lambda){
    return(pricefutures(ttm = x[-(1:2)], s0 = x[1], delta0 = x[2],
                        sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                        sigmaE = sigmaE, rho = rho, r = r, lambda = lambda))
  }

  p.futures <- t(apply(cbind(state$state, ttm), 1, pricefutures.wrapper,
                       sigmaS = coefs$sigmaS, kappa = coefs$kappa,
                       alpha = coefs$alpha, sigmaE = coefs$sigmaE,
                       rho = coefs$rho, r = coefs$r, lambda = coefs$lambda))

  dimnames(p.futures) <- dimnames(data)
  return(p.futures)
}
### <---------------------------------------------------------------------->
setMethod("fitted", signature(object = "schwartz2f.fit"), fitted.schwartz2f.fit)
### <---------------------------------------------------------------------->


### <======================================================================>
vcov.schwartz2f <- function(object, time = 1)
{

  if(length(time) != 1){
    stop("'time' must be a scalar!")
  }
  
  tmp.coef <- coef(object)

  sigma.log <- .sigma.state.schwartz2f(sigmaS = tmp.coef$sigmaS,
                                       kappa = tmp.coef$kappa,
                                       sigmaE = tmp.coef$sigmaE,
                                       rho = tmp.coef$rho,
                                       time = time)
  return(sigma.log)
}
### <---------------------------------------------------------------------->
setMethod("vcov", signature(object = "schwartz2f"), vcov.schwartz2f)
### <---------------------------------------------------------------------->


### <======================================================================>
mean.schwartz2f <- function(x, time = 1)
{

  time <- .get.data(time, type = "uv")
  
  .mean <- function(time, object){
    tmp.coef <- coef(object)
    means.log <- .mu.state.schwartz2f(x0 = log(tmp.coef$s0),
                                      delta0 = tmp.coef$delta0,
                                      mu = tmp.coef$mu,
                                      sigmaS = tmp.coef$sigmaS,
                                      kappa = tmp.coef$kappa,
                                      alpha = tmp.coef$alpha,
                                      sigmaE = tmp.coef$sigmaE,
                                      rho = tmp.coef$rho,
                                      time = time)

    sigma.log <- vcov(x, time = time)
    return(c(s.t = exp(means.log[1] + 1/2 * sigma.log[1,1]),
             delta.t = means.log[2]))

  }
  
  ## Vectorize the mean function
  means <- t(apply(matrix(time), 1, .mean, object = x))

  if(length(time) == 1){
    means <- c(s.t = means[1], delta.t = means[2])
  }
  
  return(means)

}
### <---------------------------------------------------------------------->
setMethod("mean", signature(x = "schwartz2f"), mean.schwartz2f)
### <---------------------------------------------------------------------->

### <======================================================================>
show.schwartz2f <- function(object)
{
    cat("\n----------------------------------------------------------------\n")
    cat("Schwartz97 two-factor model:\n\n")
    cat("SDE\n")
    cat("d S_t     = S_t   * (mu - delta_t)    * dt + S_t * sigmaS * dW_1\n")
    cat("d delta_t = kappa * (alpha - delta_t) * dt + sigmaE * dW_2\n")
    cat("E(dW_1 * dW_2) = rho * dt\n\n")
    tmp.coef <- coef(object)
    ## tmp.coef.formatted <- sapply(tmp.coef, function(x)sprintf("% .3E", x))
    cat("Parameters\n")
    tmp.coef.formatted <- sapply(tmp.coef, as.character)
    invisible(apply(cbind(names(tmp.coef), tmp.coef.formatted), 1,
                    function(x, max.len)cat(x[1], rep(" ", max.len - nchar(x[1])), ": ", x[2], "\n", sep = ""),
                    max.len = max(nchar(names(tmp.coef.formatted)))))
    cat("----------------------------------------------------------------\n")
}
### <---------------------------------------------------------------------->
setMethod("show", signature(object = "schwartz2f"), show.schwartz2f)
### <---------------------------------------------------------------------->


### <======================================================================>
show.schwartz2f.fit <- function(object)
{
    cat("\n----------------------------------------------------------------\n")
    cat("Fitted Schwartz97 two-factor model:\n\n")
    cat("SDE (P-dynamcis)\n")
    cat("d S_t     = S_t   * (mu - delta_t) * dt    + S_t * sigmaS * dW_1\n")
    cat("d delta_t = kappa * (alpha - delta_t) * dt + sigmaE * dW_2\n")
    cat("E(dW_1 * dW_2) = rho * dt\n\n")
    cat("SDE (Q-dynamcis)\n")
    cat("d S_t     = S_t   * (r - delta_t) * dt      + S_t * sigmaS * dW*_1\n")
    cat("d delta_t = kappa * (alphaT - delta_t) * dt + sigmaE * dW*_2\n")
    cat("alphaT = alpha - lambda/kappa\n\n")

    tmp.coef <- coef(object)
    ## tmp.coef.formatted <- sapply(tmp.coef, function(x)sprintf("% .3E", x))
    cat("Parameters\n")
    tmp.coef.formatted <- sapply(tmp.coef, as.character)
    invisible(apply(cbind(names(tmp.coef), tmp.coef.formatted), 1,
                    function(x, max.len)cat(x[1], rep(" ", max.len - nchar(x[1])), ": ", x[2], "\n", sep = ""),
                    max.len = max(nchar(names(tmp.coef.formatted)))))
    cat("----------------------------------------------------------------\n")
    cat("Optimization information")
    cat("\nConverged:          ", object@converged, sep = "")
    cat("\nFitted parameters:  ")
    fitted.params <- names(object@fitted.params)[object@fitted.params]
    cat(paste(fitted.params, collapse = ", "), "; (Number: ", length(fitted.params), ")", sep = "")
    cat("\nlog-Likelihood:     ", object@llh, sep = "")
    cat("\nNbr. of iterations: ", object@n.iter, sep = "")
    cat("\n----------------------------------------------------------------\n")
}
### <---------------------------------------------------------------------->
setMethod("show", signature(object = "schwartz2f.fit"),
          show.schwartz2f.fit)
### <---------------------------------------------------------------------->


### <======================================================================>
coef.schwartz2f <- function(object)
{
    return(list(s0 = object@s0,
                delta0 = object@delta0,
                mu = object@mu,
                sigmaS = object@sigmaS,
                kappa = object@kappaE,
                alpha = object@alpha,
                sigmaE = object@sigmaE,
                rho = object@rhoSE))
}
### <---------------------------------------------------------------------->
setMethod("coef", signature(object = "schwartz2f"),
          coef.schwartz2f)
### <---------------------------------------------------------------------->
setMethod("coefficients", signature(object = "schwartz2f"),
          coef.schwartz2f)
### <---------------------------------------------------------------------->

### <======================================================================>
coef.schwartz2f.fit <- function(object)
{
    return(list(s0 = object@s0,
                delta0 = object@delta0,
                mu = object@mu,
                sigmaS = object@sigmaS,
                kappa = object@kappaE,
                alpha = object@alpha,
                sigmaE = object@sigmaE,
                rho = object@rhoSE,
                r = object@r,
                lambda = object@lambdaE,
                alphaT = object@alphaT))
}
### <---------------------------------------------------------------------->
setMethod("coef", signature(object = "schwartz2f.fit"),
          coef.schwartz2f.fit)
### <---------------------------------------------------------------------->
setMethod("coefficients", signature(object = "schwartz2f.fit"),
          coef.schwartz2f.fit)
### <---------------------------------------------------------------------->

### <======================================================================>
plot.schwartz2f <- function(x, n = 100, time = 2, dt = 1/52)
{

  ## trajectories <- lapply(1:n, function(dummy, obj, n, t)simstate(n, t, obj),
  ##                        obj = x, n = time/dt, t = time)
  trajectories <- replicate(n, simstate(time/dt, time, x), simplify = FALSE)
  
  st <- sapply(trajectories, function(x)x[,1])
  deltat <- sapply(trajectories, function(x)x[,2])

  time.seq <- seq(dt, time, by = dt)
  means <- mean(x, time = time.seq)     # Calculate expectations

  quantiles <- c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)
  names(quantiles) <- paste(quantiles * 100, "%", sep = "")

  log.st.var <- sapply(time.seq, function(t,obj)vcov(obj, t)[1,1], obj = x)
  st.quantiles <- sapply(quantiles, qlnorm,
                         meanlog = log(means[,1]) - 1 / 2 * log.st.var,
                         sdlog = sqrt(log.st.var))

  deltat.sd <- sqrt(sapply(time.seq, function(t,obj)vcov(obj, t)[2,2], obj = x))
  deltat.quantiles <- sapply(quantiles, qnorm, mean = means[,2], sd = deltat.sd)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
  })

  par(mfrow = c(2, 1))
  par(oma = c(5,0,0,0) + 0.1)
  par(mai = c(0,1,0,0))

  ## plot spot prices
  plot(time.seq, time.seq, ylim = range(st, st.quantiles), type = "n",
       main = "", xlab = "", ylab = "S(t)", xaxt = "n")

  apply(st, 2, function(y, x)lines(x, y, col = "grey"), x = time.seq) 
  lines(time.seq, means[,1], col = "red")
  lty <- c("dotted", "dashed", "dotdash")
  for(i in 1:6)
    lines(time.seq, st.quantiles[,i], col = "red", lty = c(lty, rev(lty))[i])
  legend("topleft",
         c("Trajectories", "Mean", "99% CI", "95% CI", "90% CI"),
         col = c("grey", rep("red", 4)), lty = c("solid", "solid", lty), bg = "white")
  
  ## plot convenience yield
  plot(time.seq, time.seq, ylim = range(deltat, deltat.quantiles), type = "n",
       main = "", xlab = "time", ylab = "delta(t)")
  apply(deltat, 2, function(y, x)lines(x, y, col = "grey"), x = time.seq) 
  lines(time.seq, means[,2], col = "red")
  for(i in 1:6)
    lines(time.seq, deltat.quantiles[,i], col = "red", lty = c(lty, rev(lty))[i])
}

### <---------------------------------------------------------------------->
setMethod("plot", signature(x = "schwartz2f", y = "missing"), plot.schwartz2f)
### <---------------------------------------------------------------------->

### <======================================================================>
plot.schwartz2f.fit <- function(x, type = c("trace.pars", "state", "forward.curve", "sim"),
                                data, ttm, ...)
{
  type <- match.arg(type)
  if(type == "trace.pars"){
    trace.pars <- x@trace.pars
    if(ncol(trace.pars) > 10){
      warning("Only the first 10 columns of 'x@trace.pars' are plotted due to the limitation in plot.ts!")
      trace.pars <- trace.pars[, 1:10]
    }
    trace.pars[,"rel.tol"] <- log(trace.pars[,"rel.tol"], base = 10)
    colnames(trace.pars)[colnames(trace.pars) == "rel.tol"] <- "log10(rel.tol)"
    plot(as.ts(trace.pars), xlab = "Iteration", type = "p",
         main = "Parameter evolution",
         panel = function(x,...){points(x,...);abline(h=rev(x)[1])}, ...)
  }else if(type == "state"){
    if(missing(data) | missing(ttm)){
      stop("'data' and 'ttm' must be submitted if type == 'state'")
    }
    state <- filter.schwartz2f(data, ttm, x)$state
    col <- colorRampPalette(c("darkblue", "lightblue"))(ncol(data))
      
    oldpar <- par(no.readonly = TRUE)
    on.exit({
      par(oldpar)
    })

    par(mfrow = c(2, 1))
    par(oma = c(5,0,0,0) + 0.1)
    par(mai = c(0,1,0,0))

    plot(as.Date(rownames(data)), state[,1], type = "l",
         main = "", xlab = "", ylab = "Spot and futures price", xaxt = "n", ...)
    for(i in 1:ncol(data))
      lines(as.Date(rownames(data)), data[,i], type = "l", col = col[i])

    legend("topleft", legend = c("Spot", "1st Futures", "Last futures"),
           fill = c("black", col[1], rev(col)[1]))
    plot(as.Date(rownames(data)), state[,2], type = "l", xlab = "", ylab = "Convenience yield", ...)
    abline(h = coef(x)$alpha)
  }else if(type == "sim"){
    callNextMethod(x, ...)
  }else if(type == "forward.curve"){

    state <- filter.schwartz2f(data, ttm, x)$state
    fitted.futures <- fitted(x, data, ttm)

    x.seq <- seq(0, by = x@deltat, length = nrow(state))

    plot(x.seq, state[,1], type = "l",
         xlim = c(min(x.seq), max(x.seq) + max(ttm[nrow(ttm),], na.rm = TRUE)),
         ylim = range(cbind(state[,1], fitted.futures), na.rm = TRUE),
         ylab = "", xlab = "Time", ...)

    plot.forward.curve <- function(i, for.cur, ttm.mat, date.idx, col, ...){
      lines(date.idx[i] + c(0, ttm.mat[i,]),  for.cur[i,],
            col = col[i], lty = "dotted", ...)
    }

    sapply(1:nrow(data), plot.forward.curve, for.cur = cbind(state[,1], fitted.futures),
           ttm.mat = ttm, col = rep(rainbow(10), length = nrow(data)),
           date.idx = x.seq, ...)
    legend("topleft", legend = c("Filtered spot price", "Fitted Forward Curves"),
           lty = c("solid", "dotted"), bg = "white")

  }
}

### <---------------------------------------------------------------------->
setMethod("plot", signature(x = "schwartz2f.fit", y = "missing"), plot.schwartz2f.fit)
### <---------------------------------------------------------------------->
