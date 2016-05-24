

profile.glm <- function(fitted, which.par, alpha = 0.005, max.steps = 50,
               nsteps = 8, step.warn = 5, trace = F, ...)
{
  ## Match and test input arguments
  call <- match.call()
  name.fitted <- call$fitted # Name of original fitted object
  if(!any(class(fitted) == 'glm'))
    stop("Argument 'fitted' must be of class glm")
  if(family(fitted)$family != 'binomial')
    stop("GLM must be fitted with the binomial family")
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(round(max.steps) > round(nsteps))
  stopifnot(round(nsteps) > round(step.warn))
  stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
  max.steps <- round(max.steps)
  nsteps <- round(nsteps)
  step.warn <- round(step.warn)
  trace <- as.logical(trace)[1]
  ## Misc extraction from 'fitted'
  call.fitted <- fitted$call # Original call
  mf <- model.frame(fitted)
  X <- model.matrix(fitted)
  y <- model.response(mf)
  n <- length(y)
  if(!is.null(dim(y))) n <- n/2
  if(!is.null(model.weights(mf))) w <- model.weights(mf)
  else w <- rep(1, n)
  if(!is.null(model.offset(mf))) o <- model.offset(mf)
  else o <- rep(0, n)
  p <- length(pnames <- names(b0 <- coef(fitted)))
  na.par <- is.na(b0)
  if(missing(which.par)) which.par <- 1:p
  if(!is.numeric(which.par) && !is.character(which.par))
    stop("'which.par' must be a character or integer vector")
  if(is.numeric(which.par)) {
    if(any(abs(which.par - round(which.par)) > 0.0001)) {
    which.par <- round(which.par)
    warning("Non-integer element(s) in 'which.par' rounded off to the nearest integer") 
    }
  }
  if(is.character(which.par)) which.par <- match(which.par, pnames, 0)
  if(any(is.na(which.par)))
    stop("Invalid parameter argument(s) in 'which.par'")
  stopifnot(length(which.par) > 0)
  std.err <- sqrt(diag(vcov(fitted)))
  orig.dev <- deviance(fitted)
  lLik <- logLik(fitted)
  fam <- family(fitted)
  ## Profile limit
  lroot.max <- qnorm(1 - alpha/2)
  ## Step length
  delta <- lroot.max / nsteps
  ## Result list
  prof.list <- vector('list', length=length(which.par))
  names(prof.list) <- pnames[which.par]
  ## for each which.par move down and then up (from MLE), fit model
  ## with 'wp' fixed at ascending/descending values (using offset) and
  ## store signed likelihood root values and parameter values
  for(wp in which.par) {
    if(na.par[wp]) next ## If coef is NA skip this and go to next wp
    skip.wp <- na.par
    skip.wp[wp] <- TRUE ## Exclude coef with NA and current wp from
                        ## design matrix
    X.wp <- X[, !skip.wp, drop = FALSE] ## Design matrix without wp
    par.wp <- b0[wp] ## Reset wp coefficient to initial value
    lroot.wp <- 0 ## lroot at MLE
    wp.name <- pnames[wp] ## Name of wp parameter
    for(direction in c(-1, 1)) { ## Move down and then up
      if(trace) {
        message("\nParameter '", wp.name, "' ",
               c('down','up')[(direction + 1)/2 + 1])
        utils::flush.console()
      }
      step <- 0
      lroot <- 0
      etastart <- X[, !na.par, drop = FALSE] %*% b0[!na.par] + o
      ## Increment wp (using offset) and refit model until threshold is
      ## reached 
      while(step < max.steps && abs(lroot) <= lroot.max) {
        step <- step + 1
        bi <- b0[wp] + direction * step * delta * std.err[wp.name]
        oi <- as.vector(o + X[, wp] * bi)
        fm <- glm.fit(x = X.wp, y = y, weights = w, etastart = etastart,
                      offset = oi, family = fam, control = fitted$control)
        dev <- (fm$deviance - orig.dev)
        ## Likelihood root statistic
        lroot <- direction * sqrt(dev) 
        lroot.wp <- c(lroot.wp, lroot)
        par.wp <- c(par.wp, bi)    
        ## Update start values
        etastart <- X[, !skip.wp, drop = FALSE] %*% fm$coefficients + oi
      } ## end 'while step <= max.steps and lroot <= lroot.max'
      ## test that lroot.max is reached (e.g. max.steps not reached)
      if(abs(lroot) < lroot.max)
        warning("profile may be unreliable for '", wp.name,
                "' because lroot.max was not reached for ",
                wp, c(" down", " up")[(direction + 1)/2 + 1])
      ## test that enough steps are taken
      if(step <= step.warn)
        warning("profile may be unreliable for '", wp.name,
                "' because only ", step, "\n  steps were taken ",
                c("down", "up")[(direction + 1)/2 + 1])
    } ## end 'for direction down or up'
    par.order <- order(par.wp)
    prof.list[[wp.name]] <-
      structure(data.frame(par.wp[par.order]), names = "par.values")
    prof.list[[wp.name]]$lroot <- lroot.wp[par.order]
    if(!all(diff(lroot.wp[par.order]) > 0))
      warning("likelihood is not monotonically decreasing from maximum,\n",
              "  so profile may be unreliable for ", wp.name)
  } ## end 'for wp in which.par'
  val <- structure(prof.list, original.fit = fitted, summary =
    summary(fitted), name.originalfit = name.fitted, call.originalfit =
    call.fitted)  
  class(val) <- c('profile.glm')
  return(val)
}



plot.profile.glm <- function(x, which.par, likelihood = TRUE,
                             log = FALSE, relative = TRUE, approx = TRUE,
                             conf.int = TRUE, level = 0.95, n = 100, 
                             fig = TRUE, ylim = NULL, ...)
{  
  ## Match and test arguments
  xnames <- names(x)
  if(missing(which.par)) which.par <- seq_along(xnames)
  if(!is.numeric(which.par) && !is.character(which.par))
    stop("'which.par' must be a character or integer vector")
  if(is.numeric(which.par)) {
    if(any(abs(which.par - round(which.par)) > 0.0001)) {
    which.par <- round(which.par)
    warning("Non-integer element(s) in 'which.par' rounded off to the nearest integer")
    }
    if(any(which.par > length(x)))
       stop("One of more elements in 'which.par' exceedes number of parameters in x")
  }
  if(is.character(which.par)) which.par <- match(which.par, xnames, 0)
  if(any(is.na(which.par)))
    stop("Invalid parameter argument(s) in 'which.par'")
  stopifnot(length(which.par) > 0)
  ## Match extraction from original 'fitted' glm
  orig.fit <- attr(x, 'original.fit')
  name.originalfit <- attr(x, 'name.originalfit')
  call.originalfit <- attr(x, 'call.originalfit') 
  b0 <- coef(orig.fit)
  std.err <- sqrt(diag(vcov(orig.fit)))
  lLik <- logLik(orig.fit)
  par.null <- sapply(x, is.null)
  if(any(par.null[which.par]))
    warning(paste("Profile likelihood not shown for parameter '",
                  names(which(par.null)),
                  "' due to missing parameter estimate", sep='')) 
  ## Result list
  res.list <- c(name.originalfit, call.originalfit,
    vector('list', length=length(which.par)))
  names(res.list) <- c('name.originalfit', 'call.originalfit',
    xnames[which.par]) 
  ## Open a device if fig = T
  if(fig) {
    nrows <- ceiling(sqrt(length(which.par[!par.null[which.par]])))
    if(nrows > 3) {
      op <- par(mfrow=c(3,3), mar = c(4.5,4,2,1), oma = c(1,1,5,1))
      oask <- devAskNewPage()
    }
    else op <- par(mfrow=c(nrows, nrows), mar = c(4.5,4,2,1), oma = c(1,1,5,1))
    on.exit(par(op))
  } 
  ## For each wp calculate requisite values and plot them  
  for(wp in which.par) {
    if(par.null[wp]) next
    par.wp <- x[[wp]]$par.values
    lik.wp <- x[[wp]]$lroot
    spline.wp <- spline(par.wp, lik.wp, n = n)
    wp.name <- xnames[wp]
    if(approx) {
      approx.wp <- (par.wp - b0[wp.name])/std.err[wp.name]
      spline.approx.wp <- spline(par.wp, approx.wp, n = n)
    }
    if(conf.int) cutoff.wp <- qnorm((1-level)/2, lower.tail = F) 
    if(!likelihood) { ## lroot
      ylab <- "Profile likelihood root"
      if(conf.int) cutoff.wp <- c(-1, 1) %o% cutoff.wp
    } ## end 'if !likelihood'
    else { ## if likelihood == TRUE
      spline.wp$y <- -spline.wp$y^2/2
      if(approx) spline.approx.wp$y <- -spline.approx.wp$y^2/2
      if(conf.int) cutoff.wp <- -cutoff.wp^2/2
      if(log && relative)
        ylab <- "Rel. profile log-likelihood"
      if(!log && relative) {
        ylab <- "Rel. profile likelihood"
        spline.wp$y <- exp(spline.wp$y)
        if(approx) spline.approx.wp$y <- exp(spline.approx.wp$y)
        if(conf.int) cutoff.wp <- exp(cutoff.wp)
      }
      if(log && !relative) {
        ylab <- "Profile log-likelihood"
        spline.wp$y <- spline.wp$y + lLik
        if(approx) spline.approx.wp$y <- spline.approx.wp$y + lLik
        if(conf.int) cutoff.wp <- cutoff.wp + lLik
      }
      if(!log && !relative) {
        ylab <- "Profile likelihood"
        spline.wp$y <- exp(spline.wp$y + lLik)
        if(approx) spline.approx.wp$y <- exp(spline.approx.wp$y + lLik)
        if(conf.int) cutoff.wp <- exp(cutoff.wp + lLik)
      }
    } ## end 'if likelihood == TRUE'
    res.list[[wp.name]]$spline.vals <- spline.wp
    if(approx) res.list[[wp.name]]$spline.approx.vals <- spline.approx.wp 
    if(conf.int) res.list[[wp.name]]$cutoff <- cutoff.wp 
    if(fig) {
      plot(spline.wp$x, spline.wp$y, type='n', xlab=wp.name,
      ylab=ylab, col='black', ylim = ylim)
      if(!likelihood) points(b0[wp.name], 0, pch=18)
      if(conf.int) abline(h = cutoff.wp)
      if(approx) {
        lines(spline.approx.wp$x, spline.approx.wp$y, lty=2,
          col='steelblue')
      }
      lines(spline.wp$x, spline.wp$y)
      title(paste("\nProfile likelihood of parameter estimates from '",
      name.originalfit, "'", sep=""), outer=T)
      if(nrows > 3) {
        devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
      }
    } ## end 'if fig = T draw plot for each wp'
  } ## end 'for wp in which.par calculate values needed for plotting'
  if(fig) invisible(res.list)
  else return(res.list)
}



