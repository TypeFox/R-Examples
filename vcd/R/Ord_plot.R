# This should be revised to allow graphical parameters to be more easily passed
# for points and lines
# For now, added lwd, lty and col args for lines, with more useful defaults

Ord_plot <- function(obj, legend = TRUE, estimate = TRUE, tol = 0.1,
                     type = NULL, xlim = NULL, ylim = NULL, xlab = "Number of occurrences",
		     ylab = "Frequency ratio", main = "Ord plot", gp = gpar(cex = 0.5),
		     lwd = c(2,2), lty=c(2,1), col=c("black", "red"),
		     name = "Ord_plot", newpage = TRUE, pop = TRUE,
                     return_grob = FALSE, ...)
{
  if(is.vector(obj)) {
    obj <- table(obj)
  }
  if(is.table(obj)) {
    if(length(dim(obj)) > 1) stop ("obj must be a 1-way table")
    x <- as.vector(obj)
    count <- as.numeric(names(obj))
  } else {
    if(!(!is.null(ncol(obj)) && ncol(obj) == 2))
      stop("obj must be a 2-column matrix or data.frame")
    x <- as.vector(obj[,1])
    count <- as.vector(obj[,2])
  }

  y <- count * x/c(NA, x[-length(x)])
  fm <- lm(y ~ count)
  fmw <- lm(y ~ count, weights = sqrt(pmax(x, 1) - 1))
  fit1 <- predict(fm, data.frame(count))
  fit2 <- predict(fmw, data.frame(count))
  if(is.null(xlim)) xlim <- range(count)
  if(is.null(ylim)) ylim <- range(c(y, fit1, fit2), na.rm = TRUE)
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.04
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.04

	lwd <- rep_len(lwd, 2)  # assure length=2
	lty <- rep_len(lty, 2)
	col <- rep_len(col, 2)

  if(newpage) grid.newpage()
  pushViewport(plotViewport(xscale = xlim, yscale = ylim, default.units = "native", name = name))
  grid.points(x = count, y = y, default.units = "native", gp = gp, ...)
  grid.lines(x = count, y = fit1, default.units = "native", gp = gpar(lwd=lwd[1], lty=lty[1], col=col[1]))
  grid.lines(x = count, y = fit2, default.units = "native", gp = gpar(lwd=lwd[2], lty=lty[2], col=col[2]))
  grid.rect(gp = gpar(fill = "transparent"))
  grid.xaxis()
  grid.yaxis()
  grid.text(xlab, y = unit(-3.5, "lines"))
  grid.text(ylab, x = unit(-3, "lines"), rot = 90)
  grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))

  RVAL <- coef(fmw)
  names(RVAL) <- c("Intercept", "Slope")
  if(legend)
  {
    legend.text <- c(paste("slope =", round(RVAL[2], digits = 3)),
                     paste("intercept =", round(RVAL[1], digits = 3)))
    if(estimate) {
      ordfit <- Ord_estimate(RVAL, type = type, tol = tol)
      legend.text <- c(legend.text, "", paste("type:", ordfit$type),
        paste("estimate:", names(ordfit$estimate),"=", round(ordfit$estimate, digits = 3)))
      legend.text <- paste(legend.text, collapse = "\n")
    }
    grid.text(legend.text, min(count), ylim[2] * 0.95, default.units = "native", just = c("left", "top"))
  }

  if(pop) popViewport() else upViewport()
  if(return_grob)
      invisible(structure(RVAL, grob = grid.grab()))
  else
      invisible(RVAL)
}

Ord_estimate <- function(x, type = NULL, tol = 0.1)
{
  a <- x[1]
  b <- x[2]
  if(!is.null(type))
    type <- match.arg(type, c("poisson", "binomial", "nbinomial", "log-series"))
  else {
    if(abs(b) < tol) type <- "poisson"
    else if(b < (-1 * tol)) type <- "binomial"
    else if(a > (-1 * tol)) type <- "nbinomial"
    else if(abs(a + b) < 4*tol) type <- "log-series"
    else type <- "none"
  }

  switch(type,

  "poisson" = {
    par <- a
    names(par) <- "lambda"
    if(par < 0) warning("lambda not > 0")
  },
  "binomial" = {
    par <- b/(b - 1)
    names(par) <- "prob"
    if(abs(par - 0.5) > 0.5) warning("prob not in (0,1)")
  },
  "nbinomial" = {
    par <- 1 - b
    names(par) <- "prob"
    if(abs(par - 0.5) > 0.5) warning("prob not in (0,1)")
  },

  "log-series" = {
    par <- b
    names(par) <- "theta"
    if(par < 0) warning("theta not > 0")
  },
  "none" = {
    par <- NA
  })
  list(estimate = par, type = type)
}


