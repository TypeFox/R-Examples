plot.tlm <-
function(x, type = c("original", "transformed", "diagnosis"), observed = FALSE, xname = "x", yname = "y", level = 0.95, ...)
 {
  if (!inherits(x, "tlm"))
     stop("argument 'x' must be of class 'tlm'")
     
  if(any(is.na(coef(x$model))))
     stop("plot is not available for models with any missing estimated coefficient")
  
  ### 'type' control:
  type <- match.arg(type)
  
  ### 'xname' control:
  if (!is.null(xname) && (!inherits(xname, "character") || is.na(xname) || length(xname) != 1))
    stop("the name for the explanatory variable X 'xname' must be a character")
  
  ### 'yname' control:
  if (!is.null(yname) && (!inherits(yname, "character") || is.na(yname) || length(yname) != 1))
    stop("the name for the response variable Y 'yname' must be a character")
  
  ### 'level' control:
  if (!inherits(level, "numeric") || level <= 0 || level >= 1 || length(level) != 1)
    stop("'level' must be a number in (0, 1)")
      
  if (type == "diagnosis")
   {
   	par(mfrow = c(2, 2), ...)
    # Diagnosis plot for the fitted model (transformed space):
    plot(x$model, ...)
    } else {
    mod <- x$model
    family <- family(mod)$family
    mf <- model.frame(mod)
    mt <- attr(mf, "terms")
    Xclass <- attr(mt, "dataClasses")[2]
    if (missing(xname)) xlabel <- names(Xclass) else xlabel <- xname
    if (missing(yname)) yname <- names(mf)[1]
    if (Xclass == "factor")
     {
      if (observed)
       warning("the observations are not shown in the plot if the explanatory variable is categorical")
      dat <- MY(x, space = type, level = level)
      M <- dat$M
      ymin <- min(M[, -1])
      ymax <- max(M[, -1])
      nlevels <- nrow(M)
      if (type == "original")
       {
       	ylabelpre <- switch(dat$ymeasure,
       	                   "geometric mean" = "Geometric mean",
                           "mean" = "Mean",
                           "median" = "Median",
                           "probability" = "Probability")
        ylabel <- paste(ylabelpre, "of", yname)
        } else {
        ylabel <- switch(family,
                         "gaussian" = paste("Mean of", yname),
                         "binomial" = paste("Log(Odds of ", yname, ")", sep = ""),
                         "poisson" = paste("Log(Mean of ", yname, ")", sep = ""))
        if (x$ypow == 0)
         ylabel <- paste("Mean of log(", yname, ")", sep = "")
        if (x$ypow != 1 & x$ypow != 0)
         ylabel <- substitute("Mean of " * ynam * phantom("")^power, list(ynam = yname, power = attr(x, "ypowlabel")))
       }
      delta <- 0.2
      plot(c(1 - delta, nlevels + delta), xlim = c(1 - delta, nlevels + delta), ylim = c(ymin, ymax), type = "n", xaxt = "n", xlab = xlabel, ylab = ylabel, ...)
      axis(1, at = 1:nlevels, labels = levels(mf[, 2]))
      lines(1:nlevels, M[, 2], lty = 2, col = "black")
      points(1:nlevels, M[, 2], pch = 19, ...)
      segments(1:nlevels, M[, 3], 1:nlevels, M[, 4], lwd = 1.5, ...)
      } else {
      dat <- MY(x, npoints = 500, space = type, level = level)
      M <- dat$M
      if (type == "original")
       {
       	ylabelpre <- switch(dat$ymeasure,
       	                   "geometric mean" = "Geometric mean",
                           "mean" = "Mean",
                           "median" = "Median",
                           "probability" = "Probability")
        ylabel <- paste(ylabelpre, "of", yname)
        } else {
        ylabel <- switch(family,
                         "gaussian" = paste("Mean of", yname),
                         "binomial" = paste("Log(Odds of ", yname, ")", sep = ""),
                         "poisson" = paste("Log(Mean of ", yname, ")", sep = ""))   
        if (x$xpow == 0)
         xlabel <- paste("Log(", xname, ")", sep = "")
        if (x$xpow != 1 & x$xpow != 0)
         xlabel <- substitute(xnam * phantom("")^power, list(xnam = xname, power = attr(x, "xpowlabel")))
        if (x$ypow == 0)
         ylabel <- paste("Mean of log(", yname, ")", sep = "")
        if (x$ypow != 1 & x$ypow != 0)
         ylabel <- substitute("Mean of " * ynam * phantom("")^power, list(ynam = yname, power = attr(x, "ypowlabel")))
       }
      if (observed)
       {
       	if (family != "gaussian")
       	 warning("the observations are not shown in the plot for models different than the linear regression model (i.e., family 'gaussian'")
       	if (family == "gaussian")
       	 {
       	  # Plot with observations:
          yobs <- model.response(mf)
          xobs <- model.frame(mod)[, 2]
          if (type == "original") 
           {
       	    if (x$ypow == 0) yobs <- exp(yobs) else yobs <- yobs^(1 / x$ypow)
       	    if (x$xpow == 0) xobs <- exp(xobs) else xobs <- xobs^(1 / x$xpow)
           }
          ymin <- min(M[, -1], yobs)
          ymax <- max(M[, -1], yobs)
          x <- xobs
          y <- yobs
          #plot(x, y, type = "p", col = "gray", pch = 19, cex = 0.6, ylim = c(ymin, ymax), xlab = xlabel, ylab = ylabel, ...)
          plot(x, y, type = "p", col = "gray", ylim = c(ymin, ymax), xlab = xlabel, ylab = ylabel, ...)
          lines(M[, 1], M[, 2], ...)
          lines(M[, 1], M[, 3], lty = 2, ...)
          lines(M[, 1], M[, 4], lty = 2, ...)
         }
        } else {
        # Plot with no observations:
         ymin <- min(M[, -1])
         ymax <- max(M[, -1])
         plot(M[, 1], M[, 2], type = "l", ylim = c(ymin, ymax), xlab = xlabel, ylab = ylabel, ...)
         lines(M[, 1], M[, 3], lty = 2, ...)
         lines(M[, 1], M[, 4], lty = 2, ...)
        }
     }
   }
 }
