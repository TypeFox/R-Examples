which.plots <- function(x, which, ask, model.name, nx, ...)
{
  args <- list(...)
  nw <- length(which)
  if((!ask || nx > 1L) && nw > 1L)
    par(mfrow = n2mfrow(nw))
  else
    ask <- TRUE  
  if(nx > 1L)
    ask <- TRUE
  if(nw > 1L && ask && nx > 1L) {
    par(oma = c(1, 1, 5, 1))
    if(is.null(args$main))
      main <- paste("Diagnostic plots for model", model.name)
    else
      main <- args$main
  }
  residuals <- x$residuals
  if(is.matrix(residuals))
    residuals <- residuals[,1L]
  residuals <- as.numeric(residuals)
  fitvalues <- x$fitted.values
  if(is.matrix(fitvalues))
    fitvalues <- fitvalues[,1L]
  fitvalues <- as.numeric(fitvalues)
  k <- 0L
  ok <- FALSE
  for(ww in which) {
    k <- k + 1L
    if(ww == "hist-resid" && length(residuals) > 0L) {
      args2 <- args
      dens <- density(residuals)
      hst <- hist(residuals, plot = FALSE)
      args2$ylim <- c(0, max(c(hst$density, dens$y)))
      args2$xlab <- "Residuals"
      args2$ylab <- "Density"
      args2$main <- "Histogramm and density"
      args2$freq <- FALSE
      args2$x <- residuals
      args2 <- delete.args(graphics::hist.default, args2)
      args2$cex <- args$cex
      do.call(graphics::hist.default, args2)
      graphics::lines(dens)
      box()
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "qq-resid" && length(residuals) > 0L) {
      args2 <- args
      r2 <- (residuals - mean(residuals))/sd(residuals)
      args2$y <- r2
      args2$ylab <- "Standardized residuals"
      args2$xlab <- "Theoretical quantiles"
      args2$main <- "Normal Q-Q Plot"
      args2$ylim <- NULL
      args2$xlim <- NULL
      args2 <- delete.args("qqnorm.default", args2, package = "stats")
      args2$cex <- args$cex
      do.call(stats::qqnorm, args2)
      stats::qqline(r2)
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "scatter-resid" && length(residuals) > 0L  && length(fitvalues) > 0L) {
      args2 <- args
      args2$y <- residuals
      args2$x <- fitvalues
      args2 <- delete.args(stats::scatter.smooth, args2)
      args2$xlab <- "Fitted values"
      args2$ylab <- "Residuals"
      args2$main <- "Fitted values vs. residuals"
      args2$cex <- args$cex
      do.call(stats::scatter.smooth, args2)
      abline(h = 0, lty = 2)
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "scale-resid" && length(residuals) > 0L  && length(fitvalues) > 0L) {
      args2 <- args
      args2$y <- sqrt(abs((residuals - mean(residuals))/sd(residuals)  ))
      args2$x <- fitvalues
      args2 <- delete.args(stats::scatter.smooth, args2)
      args2$xlab <- "Fitted values"
      args2$ylab <- expression(sqrt(abs("Standardized residuals")))
      args2$main <- "Scale-location" 
      args2$cex <- args$cex
      do.call(stats::scatter.smooth, args2)
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "scale-samples") {
      if(!is.null(attr(x$variance, "sample"))) {
        args2 <- args
        args2$x <- attr(x$variance, "sample")
        args2$selected <- "scale"
        args2$var <- TRUE
        if(is.null(args2$main))
          args2$main <- "Variance sampling path"
        do.call("plotsamples", args2)	
      }
    }
  }
  if(nw > 1L && ask && nx > 1L && ok)
    mtext(main, side = 3L, line = 2L, outer = TRUE, font = 2, cex = 1)
  if(k == 1L && nx > 1L)
    par(ask = ask) 
  if(!ok)
    warning("there is nothing to plot!")

  return(invisible(NULL))
}

