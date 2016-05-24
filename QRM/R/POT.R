## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


findthreshold <- function(data, ne)
{
    if(is.timeSeries(data)) data <- as.vector(series(data))
    if(!is.vector(data))
        stop("\ndata input to findthreshold() must be a vector or timeSeries with only one data column.\n")
    if(all(length(data) < ne)) stop("\ndata length less than ne (number of exceedances.\n")
    data <- rev(sort(as.numeric(data)))
    thresholds <- unique(data)
    indices <- match(data[ne], thresholds)
    indices <- pmin(indices + 1., length(thresholds))
    thresholds[indices]
}

## TODO: return "data" are actually the exceedances, not the input data!!!
##       => bad naming here (use 'x' for input)
fit.GPD <- function(data, threshold = NA, nextremes = NA, type = c("ml", "pwm"),
                    information = c("observed", "expected"),
                    optfunc = c("optim", "nlminb"), verbose=TRUE, ...)
{
  type <- match.arg(type)
  optfunc <- match.arg(optfunc)
  information <- match.arg(information)
  if(is.na(nextremes) & is.na(threshold))
    stop("\nEnter either a threshold or the number of upper extremes.\n")
  if(!is.na(nextremes) & !is.na(threshold))
    stop("\nEnter either a threshold or the number of upper extremes.\n")
  if(is.timeSeries(data)) data <- series(data)
  data <- as.numeric(data)
  n <- length(data)
  if(!is.na(nextremes)){
    threshold <- findthreshold(data, nextremes)
  }
  exceedances <- data[data > threshold]
  excess <- exceedances - threshold
  Nu <- length(excess)

  ## probability weighted moments (pwm)
  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0.0
  pvec <- ((1.:Nu) + delta)/(Nu + delta)
  a1 <- mean(sort(excess) * (1. - pvec))
  xi <- 2. - a0/(a0 - 2. * a1) # estimated xi
  beta <- (2. * a0 * a1)/(a0 - 2. * a1) # estimated beta
  par.ests <- c(xi, beta) # initial estimates (improved by type = "ml")

  if(type == "ml"){ # maximum likelihood (based on pwm initial values)
    negloglik <- function(theta, ydata) -sum(dGPD(ydata, theta[1], abs(theta[2]), log = TRUE))
    deriv <- function(theta, ydata){
      xi <- theta[1]
      beta <- theta[2]
      term1 <- sum(ydata / (beta + xi * ydata))
      term2 <- sum(log(1 + xi * ydata / beta))
      d1 <- -term2 * xi^(-2) + (1 + 1 / xi) * term1
      d2 <- (length(ydata) - (xi + 1) * term1) / beta
      c(d1, d2)
    }
    if(optfunc == "optim"){
      fit <- optim(par.ests, fn = negloglik, gr = deriv, ydata = excess, ...)
    }
    if(optfunc == "nlminb"){
      fit <- nlminb(start = par.ests, objective = negloglik, gradient = deriv, ydata = excess, ...)
    }
    par.ests <- fit$par
    par.ests[2] <- abs(par.ests[2])
    ifelse(fit$convergence == 0, converged <- TRUE, converged <- FALSE)
    ll.max <- -negloglik(fit$par, ydata = excess)
    if(information == "observed"){
      fisher <- hessian(negloglik, fit$par, ydata = excess)
      varcov <- solve(fisher)
    }
    if(information == "expected"){
      one <- (1. + par.ests[1.])^2. / Nu
      two <- (2. * (1. + par.ests[1.]) * par.ests[2.]^2.) / Nu
      cov <-  - ((1. + par.ests[1.]) * par.ests[2.]) / Nu
      varcov <- matrix(c(one, cov, cov, two), 2.)
    }
  }

  if(type == "pwm"){ # probability weighted moments
    denom <- Nu * (1. - 2. * xi) * (3. - 2. * xi)
    if(xi > 0.5){
      denom <- NA
      if(verbose) warning("\nAsymptotic standard errors not available for PWM Type when xi > 0.5.\n")
    }
    one <- (1. - xi) * (1. - xi + 2. * xi^2.) * (2. - xi)^2.
    two <- (7. - 18. * xi + 11. * xi^2. - 2. * xi^3.) * beta^2.
    cov <- beta * (2. - xi) * (2. - 6. * xi + 7. * xi^2. - 2. *xi^3.)
    varcov <- matrix(c(one, cov, cov, two), 2.)/denom
    information <- "expected"
    converged <- NA
    ll.max <- NA
  }

  par.ses <- sqrt(diag(varcov))
  p.less.thresh <- 1. - Nu / n
  out <- list(n = length(data), data = exceedances, threshold = threshold, # TODO: data are exceedances! not original data! confusing...
              p.less.thresh = p.less.thresh, n.exceed = Nu, type = type,
              par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              information = information, converged = converged, ll.max = ll.max)
  names(out$par.ests) <- c("xi", "beta")
  names(out$par.ses) <- c("xi", "beta")
  out
}

##' @title Plot Estimated Tail Probabilities
##' @param object object as returned by fit.GPD()
##' @param ppoints.gpd (p)points in (0,1) for evaluating GPD tail estimate
##' @param xlab x axis label
##' @param ylab y axis label
##' @param ... additional arguments passed to plot()
##' @return invisible()
##' @author Marius Hofert
plotTail <- function(object, ppoints.gpd = ppoints(256),
                     main = "Estimated tail probabilities",
                     xlab = "Exceedances x", ylab = expression(1-hat(F)[n](x)), ...)
{
    ## checks
    stopifnot(0 < ppoints.gpd, ppoints.gpd < 1)

    ## extract results from object
    exceedance <- object$data # exceedances
    u <- object$threshold # threshold
    xi.hat <- object$par.ests["xi"] # xi estimate
    beta.hat <- object$par.ests["beta"] # beta estimate

    ## compute empirical tail estimator
    ## Note: Let X be a loss and Y an excess over u (Y = X - u for X > u)
    ##       Then EKM (1997, (6.42)) says \bar{F}(u+y) = \bar{F}(u) * \bar{F}_u(y).
    ##       With estimators: \hat{\bar{F}}_n(u+y) = (Nu/n) * \hat{\bar{F}}_{u,n}(y)
    ##       If y are *sorted* excesses then u+y = sorted_exceedance and hence
    ##       \hat{\bar{F}}_n(sorted_exceedance) =  (Nu/n) * \hat{\bar{F}}_{u,n}(sorted_excess)
    ##                                          =  (Nu/n) * (Nu:1)/Nu
    ##                                          ~= (Nu/n) * rev(ppoints(Nu))
    Nu <- length(exceedance) # number of exceedances
    X. <- sort(exceedance) # sorted exceedances (sorted X which are > u)
    Nu.n <- 1 - object$p.less.thresh # empirical probability of (exceedances) x > u = Nu/n
    Fnbar.hat.X. <- Nu.n * rev(ppoints(Nu)) # \hat{\bar{F}}_n(X.) where \hat{F}_n is based on *all* data (standard edf)

    ## compute GPD tail estimator
    ## Note: As before, \hat{\bar{F}}_n(u+y) ~= (Nu/n) * rev(ppoints(Nu))
    ##       By the GPD approxi, the left-hand side equals
    ##       (Nu/n) * (1-F_{GPD(xi.hat,beta.hat)}(y))
    ##       Solving for y leads y ~= F_{GPD(xi.hat,beta.hat)}^{-1}(1-rev(ppoints(Nu)))
    ##                              = F_{GPD(xi.hat,beta.hat)}^{-1}(ppoints(Nu))
    ##       => sorted_exceedance = u+y = u + F_{GPD(xi.hat,beta.hat)}^{-1}(ppoints(Nu))
    x.gpd.est <- u + qGPD(ppoints.gpd, xi.hat, beta.hat) # 'sorted_exceedance'
    y.gpd.est <- Nu.n * rev(ppoints.gpd) # as before

    ## plot
    xlim. <- range(X., x.gpd.est)
    ## xlim. = range(object$data, x.gpd.est)
    ##       = range(object$data, u + qGPD(ppoints.gpd, xi.hat, beta.hat)) # => important for showRM()
    ylim. <- range(Fnbar.hat.X., y.gpd.est) # y range
    plot(X., Fnbar.hat.X., xlim = xlim., ylim = ylim., # make sure everything is visible
         log = "xy", main = main, xlab = xlab, ylab = ylab, ...) # plot emp. tail estimator
    lines(x.gpd.est, y.gpd.est) # add GPD tail estimator

    ## return
    invisible()
}

##' @title Plot Estimated Tail Probabilities with Risk Measure Estimates and CIs
##' @param object object as returned by fit.GPD()
##' @param alpha RM 'confidence' level (e.g., 0.999)
##' @param RM risk measure (character string)
##' @param like.num number of likelihood evaluation points for CIs
##' @param ppoints.gpd (p)points in (0,1) for evaluating GPD tail estimate
##' @param xlab x axis label
##' @param ylab y axis label
##' @param legend.pos position as accepted by legend() or NULL (no legend)
##' @param pre.0.4.9 logical indicating whether 'old return behavior' applies
##' @param ... additional arguments passed to optim()
##' @return invisible() or (if pre.0.4.9) confidence intervals
##' @author Marius Hofert
##' TODO: we should use optimize(), improve on nLL(), and use ... for plot!
showRM <- function(object, alpha, RM = c("VaR", "ES"),
                   like.num = 64, ppoints.gpd = ppoints(256),
                   xlab = "Exceedances x", ylab = expression(1-hat(F)[n](x)),
                   legend.pos = "topright", pre.0.4.9=FALSE, ...)
{
    ## checks
    stopifnot(0 < alpha, alpha < 1, like.num > 0)

    ## extract results from object
    exceedance <- object$data # exceedances
    u <- object$threshold # threshold
    xi.hat <- object$par.ests["xi"] # xi estimate
    beta.hat <- object$par.ests["beta"] # beta estimate
    Nu.n <- 1-object$p.less.thresh # empirical probability of (exceedances) x > u = Nu/n

    ## risk measure estimates
    a <- (1 - alpha) / Nu.n
    VaR <- u + beta.hat * (a^(-xi.hat) - 1) / xi.hat
    ES <- VaR / (1 - xi.hat) + (beta.hat - xi.hat * u) / (1 - xi.hat)
    RM <- match.arg(RM)
    RM.hat <- switch(RM,
                     "VaR" = VaR,
                     "ES" = ES,
                     stop("\nwrong argument 'RM'.\n")) # risk measure estimate

    ## x values for CI evaluation
    start <- switch(RM,
                    "VaR" = u,
                    "ES" = VaR,
                    stop("\nwrong argument 'RM'.\n"))
    xlim.plotTail <- range(exceedance, u + qGPD(ppoints.gpd, xi=xi.hat, beta=beta.hat)) # x range (as for plotTail()!)
    x <- exp(seq(log(start), log(xlim.plotTail[2]), length=like.num)) # x values where the evaluate likelihood (based confidence intervals)

    ## computing CIs
    nLL <- function(xi., x) { # xi. = running xi; x = x-values (risk measure value on x-axis) where the evaluate CIs
        beta <- switch(RM,
                     "VaR" = xi. * (x - u) / (a^(-xi.) - 1),
                     "ES" = ((1 - xi.) * (x - u)) / (((a^( - xi.) - 1) / xi.) + 1),
                     stop("\nwrong argument 'RM'.\n"))
        if(beta <= 0) 1e17 # dGPD() not defined there since log(beta)... but optim() needs finite initial value
        else -sum(dGPD(exceedance-u, xi., beta, log=TRUE)) # -log-likelihood
    }
    ## 1d-optimization in xi for fixed x-value (= risk measure value; each of the above 'x').
    ## The optimal *value* gives the likelihood value at the optimal xi,
    ## i.e., the y-value of the likelihood curve at the particular x (= risk measure value)
    ## (xi.hat = initial xi value)
    opt <- unlist(lapply(x, function(x.) -optim(xi.hat, fn=nLL, x=x., ...)$value))
    ## TODO more docu here
    crit <- object$ll.max - qchisq(0.999, df=1) / 2
    x <- x[opt > crit]
    opt <- opt[opt > crit]
    bothFine <- !is.na(x) & !is.na(opt)
    smth <- spline(x[bothFine], opt[bothFine], n = 256)

    ## plot
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    par(mar=c(5.1, 4.1+0.5, 4.1, 2.1+2.2)) # enlarge plot region (a bit to the left and more to the right)
    plotTail(object, ppoints.gpd=ppoints.gpd, xlab=xlab, ylab=ylab) # estimated tail probabilities (empirical and GPD)
    abline(v = RM.hat, lty=2) # vertical line indicating RM estimate
    par(new = TRUE)
    plot(x, opt, type = "n", xlab = "", ylab = "", axes = FALSE,
         xlim = xlim.plotTail, log = "x") # does *not* plot, only sets up the coordinate system (as for plotTail()!)
    lines(smth, lty=3, lwd=1.5) # RM CIs
    axis(4, at = object$ll.max - qchisq(c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999), df=1.)/2.,
         labels = c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999)) # y axis on the right
    mtext(text="Confidence intervals", side=4, line=3) # label
    if(!is.null(legend.pos)) {
        legend(legend.pos, inset=0.01, cex=0.8, lty=1:3, lwd=c(1, 1, 1.5), col=1, bty="n",
               legend=as.expression(
                      c(substitute("GPD("*xi.*","~beta.*") tail estimate",
                                   list(xi.=round(xi.hat,2), beta.=round(beta.hat,2))),
                        substitute(RM.[a]~"estimate (="~RM..*")",
                                   list(a=alpha, RM.=RM, RM..=round(RM.hat,2))),
                        substitute(RM.[a]~"CIs", list(a=alpha, RM.=RM)) )))
    }

    ci <- smth$x[smth$y > object$ll.max - qchisq(0.95, df=1)/2]
    ## return
    if(pre.0.4.9) c("Lower CI" = min(ci), "Estimate" = RM.hat, "Upper CI" = max(ci)) else invisible()
}

## ME plot
MEplot <- function(data, omit = 3., main = "Mean-Excess Plot", xlab = "Threshold", ylab = "Mean Excess", ...)
{
  if(is.timeSeries(data)) data <- series(data)
  data <- as.numeric(data)
  n <- length(data)
  myrank <- function(x, na.last = TRUE){
    ranks <- sort.list(sort.list(x, na.last = na.last))
    if(is.na(na.last)) x <- x[!is.na(x)]
    for(i in unique(x[duplicated(x)])){
      which <- x == i & !is.na(x)
      ranks[which] <- max(ranks[which])
    }
    ranks
  }
  data <- sort(data)
  n.excess <- unique(floor(length(data) - myrank(data)))
  points <- unique(data)
  nl <- length(points)
  n.excess <- n.excess[ - nl]
  points <- points[ - nl]
  excess <- cumsum(rev(data))[n.excess] - n.excess * points
  y <- excess/n.excess
  plot(points[1.:(nl - omit)], y[1.:(nl - omit)], main = main, xlab = xlab, ylab = ylab, ...)
  return(invisible(list(x = points[1.:(nl - omit)], y = y[1.:(nl - omit)])))
}

## Risk Measures
RiskMeasures <- function(out, p)
{
  u <- out$threshold
  par.ests <- out$par.ests
  xihat <- par.ests[names(par.ests) == "xi"]
  betahat <- par.ests[names(par.ests) == "beta"]
  p.less.thresh <- out$p.less.thresh
  lambda <- 1. / (1. - p.less.thresh)
  quant <- function(pp, xi, beta, u, lambda){
    a <- lambda * (1. - pp)
    u + (beta * (a^( - xi) - 1.))/xi
  }
  short <- function(pp, xi, beta, u, lambda){
    a <- lambda * (1. - pp)
    q <- u + (beta * (a^( - xi) - 1.))/xi
    (q * (1. + (beta - xi * u)/q))/(1. - xi)
  }
  q <- quant(p, xihat, betahat, u, lambda)
  es <- short(p, xihat, betahat, u, lambda)
  cbind(p, quantile = q, sfall = es)
}

## GPD shape varies with threshold or number of extremes.
xiplot <- function(data, models = 30., start = 15., end = 500., reverse = TRUE,
                   ci = 0.95, auto.scale = TRUE, labels = TRUE, table = FALSE, ...)
{
  if(is.timeSeries(data)) data <- series(data)
  data <- as.numeric(data)
  qq <- 0.
  if(ci) qq <- qnorm(1. - (1. - ci)/2.)
  x <- trunc(seq(from = min(end, length(data)), to = start, length = models))
  gpd.dummy <- function(nex, data){
    out <- fit.GPD(data = data, nextremes = nex, information = "expected")
    c(out$threshold, out$par.ests[1.], out$par.ses[1.])
  }
  mat <- apply(as.matrix(x), 1., gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("threshold", "shape", "se", "exceedances"), NULL)
  thresh <- mat[1,  ]
  y <- mat[2,  ]
  yrange <- range(y)
  if(ci){
    u <- y + mat[3.,  ] * qq
    l <- y - mat[3.,  ] * qq
    yrange <- range(y, u, l)
  }
  index <- x
  if(reverse) index <-  - x
  if(auto.scale){
    plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  } else {
    plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  }
  axis(1, at = index, labels = paste(x), tick = FALSE)
  axis(2)
  axis(3, at = index, labels = paste(format(signif(thresh, 3.))), tick = FALSE)
  box()
  if(ci){
    lines(index, u, lty = 2., col = 2.)
    lines(index, l, lty = 2., col = 2.)
  }
  if(labels){
    labely <- "Shape (xi)"
    if(ci) labely <- paste(labely, " (CI, p = ", ci, ")", sep = "")
    title(xlab = "Exceedances", ylab = labely)
    mtext("Threshold", side = 3., line = 3.)
  }
  if(table) print(mat)
  return(invisible(list(x = index, y = y, upper = u, lower = l)))
}

##' @title Hill estimator for the tail index alpha (or 1/alpha) when
##'        1-F(x) ~ x^{-alpha}*L(x), alpha > 0
##' @param data data
##' @param k number of upper order statistics
##' @param tail.index logical indicating whether the estimator of the tail index
##'        alpha is returned
##' @return Hill estimator for alpha (the default) or 1/alpha
##' @author Marius Hofert
##' @note - See EKM (1997, p. 190, Eq. (4.12))
##'       - vectorized in k
hill <- function(data, k, tail.index=TRUE)
{
    stopifnot(k >= 2, length(data) >= max(k), data > 0)
    lxs <- sort(log(data), decreasing=TRUE)
    res <- vapply(k, function(k.) mean(lxs[seq_len(k.-1)]) - lxs[k.], NA_real_)
    if(tail.index) 1/res else res
}

## Hill Plot (not using hill())
hillPlot <- function (data, option = c("alpha", "xi", "quantile"), start = 15,
    end = NA, reverse = FALSE, p = NA, ci = 0.95, auto.scale = TRUE, labels = TRUE, ...)
{
  if(is.timeSeries(data)) data <- as.vector(series(data))
  data <- as.numeric(data)
  ordered <- rev(sort(data))
  ordered <- ordered[ordered > 0]
  n <- length(ordered)
  option <- match.arg(option)
  if((option == "quantile") && (is.na(p)))
    stop("\nInput a value for the probability p.\n")
  if ((option == "quantile") && (p < 1 - start/n)){
    cat("Graph may look strange !! \n\n")
    cat(paste("Suggestion 1: Increase `p' above", format(signif(1 - start/n, 5)), "\n"))
    cat(paste("Suggestion 2: Increase `start' above ", ceiling(length(data) * (1 - p)), "\n"))
  }
  k <- 1:n
  loggs <- logb(ordered)
  avesumlog <- cumsum(loggs)/(1:n)
  xihat <- c(NA, (avesumlog - loggs)[2:n])
  alphahat <- 1/xihat
  y <- switch(option, alpha = alphahat, xi = xihat, quantile =
              ordered * ((n * (1 - p))/k)^(-1/alphahat))
  ses <- y/sqrt(k)
  if(is.na(end)) end <- n
  x <- trunc(seq(from = min(end, length(data)), to = start))
  y <- y[x]
  ylabel <- option
  yrange <- range(y)
  if(ci && (option != "quantile")){
    qq <- qnorm(1 - (1 - ci)/2)
    u <- y + ses[x] * qq
    l <- y - ses[x] * qq
    ylabel <- paste(ylabel, " (CI, p =", ci, ")", sep = "")
    yrange <- range(u, l)
  }
  if(option == "quantile") ylabel <- paste("Quantile, p =", p)
  index <- x
  if(reverse) index <- -x
  if(auto.scale){
    plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  } else {
    plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  }
  axis(1, at = index, labels = paste(x), tick = FALSE)
  axis(2)
  threshold <- findthreshold(data, x)
  axis(3, at = index, labels = paste(format(signif(threshold, 3))), tick = FALSE)
  box()
  if(ci && (option != "quantile")){
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if(labels){
    title(xlab = "Order Statistics", ylab = ylabel)
    mtext("Threshold", side = 3, line = 3)
  }
  return(invisible(list(x = index, y = y)))
}

## QQ-Plot for GPD
plotFittedGPDvsEmpiricalExcesses <- function(data, threshold = NA, nextremes = NA)
{
  if(is.na(nextremes) & is.na(threshold))
    stop("\nEnter either a threshold or the number of upper extremes.\n")
  if(is.timeSeries(data)) data <- series(data)
  mod <- fit.GPD(data, threshold, nextremes)
  if(!is.na(nextremes)) threshold <- findthreshold(as.vector(data), nextremes)
  pECDF <- edf(mod$data)
  maxVal <- as.numeric(max(data))
  quantVector <- seq(threshold, maxVal, 0.25)
  pG <- pGPD(quantVector-threshold, mod$par.ests["xi"], mod$par.ests["beta"])
  plot(quantVector, pG, type = "l", log = "x", xlab = "x (log scale)", ylab = "Fu(x-u)")
  points(mod$data, pECDF, pch = 19, col = "blue")
}



