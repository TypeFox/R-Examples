plot.betareg <- function(x, which = 1:4,
  caption = c("Residuals vs indices of obs.", "Cook's distance plot",
    "Generalized leverage vs predicted values", "Residuals vs linear predictor", 
    "Half-normal plot of residuals", "Predicted vs observed values"),
    sub.caption = paste(deparse(x$call), collapse = "\n"), main = "", 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., type = "sweighted2", nsim = 100, level = 0.9)
{
  if(!is.numeric(which) || any(which < 1) || any(which > 6)) 
    stop("`which' must be in 1:6")
    
  types <- c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2")
  Types <- c("Pearson residuals", "Deviance residuals", "Raw response residuals",
    "Weighted residuals", "Standardized weighted residuals", "Standardized weighted residuals 2")
  type <- match.arg(type, types)
  Type <- Types[type == types]

  res <- residuals(x, type = type)
  n <- length(res)
  k <- length(x$coefficients$mean)
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  Main <- rep("", 6)
  Main[which] <- rep(main, length.out = sum(show))
  one.fig <- prod(par("mfcol")) == 1
  if(ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if(show[1]) {
    plot(1:n, res, xlab = "Obs. number", ylab = Type, main = Main[1], ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[1], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[2]) {
    plot(1:n, cooks.distance(x),
      xlab = "Obs. number", ylab = "Cook's distance", type = "h", main = Main[2])
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[2], 3, 0.25)
  }
  if(show[3]) {
    plot(fitted(x), gleverage(x),
      xlab = "Predicted values", ylab = "Generalized leverage", main = Main[3], ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25)
  }
  if(show[4]) {
    plot(predict(x, type = "link"), res,
      xlab = "Linear predictor", ylab = Type, main = Main[4], ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[4], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[5]) {
    hn <- halfnormal.betareg(x, nsim = nsim, level = level, type = type)
    plot(hn[,1], hn[,2], ylim = range(hn[,-1]), main = Main[5],
      xlab = "Normal quantiles", ylab = paste(Type, "(absolute values)"), ...)
    lines(hn[,1], hn[,3],lty = 2)
    lines(hn[,1], hn[,4],lty = 1)
    lines(hn[,1], hn[,5],lty = 1)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[5], 3, 0.25)
  }
  if(show[6]) {
    y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    plot(y, fitted(x),
      xlab = "Observed values", ylab = "Predicted values", main = Main[6], ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[6], 3, 0.25)
    abline(0, 1, lty = 2, col = "gray")
  }

  if(!one.fig && par("oma")[3] >= 1) mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}

halfnormal.betareg <- function(model, nsim = 100, level = 0.90, type = "sweighted2")
{
  ## extract response y and regressors X
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model, model = "mean") else model$x$mean
  z <- if(is.null(model$x)) model.matrix(model, model = "precision") else model$x$precision
  if(is.null(model$offset$mean)) model$offset$mean <- rep(0, NROW(x))
  if(is.null(model$offset$precision)) model$offset$precision <- rep(0, NROW(z))
  wts <- weights(model)

  n <- NROW(x)
  alpha <- (1 - level)/2
  mu <- fitted(model)
  phi <- predict(model, type = "precision")
  res <- residuals(model, type = type)

  e <- matrix(0, n, nsim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:nsim) {
    ysim <- rbeta(n, mu * phi, (1 - mu) * phi)
    ctrl <- model$control
    ctrl$hessian <- FALSE
    ctrl$start <- model$coefficients
    fit <- suppressWarnings(betareg.fit(x, ysim, z, weights = wts, offset = model$offset,
      link = model$link$mean$name, link.phi = model$link$precision$name,
        control = ctrl))
    fit$y <- ysim
    fit$x <- list(mean = x, precision = z)
    class(fit) <- "betareg"
    e[,i] <- sort(abs(residuals(fit, type = type)))
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo, alpha)
    e2[i] <- quantile(eo, 1 - alpha)
  }
  
  e0 <- apply(e, 1, median)
  qq <- qnorm((n + 1:n + 0.5)/(2 * n + 1.125))
  
  cbind(qq, sort(abs(res)), e0, e1, e2)  
}
