plot.hyperblm <- function(x, breaks = "FD",
                          plotTitles = c("Residuals vs Fitted Values",
                                         "Histogram of residuals",
                                         "Log-Histogram of residuals",
                                         "Q-Q Plot"),
                          ...) {

  if (! "hyperblm" %in% class(x))
    stop("Object must be of class hyperblm")


  if(is.null(x$coef))
    stop("Object has NULL output")

  par(mar = c(6, 4, 4, 2) + 0.1)

  ## Caused problems when first argument changed to x, not needed
  ##x <- x$xMatrix
  ##y <- as.numeric(x$yVec)
  fitted.values <- x$fitted.values
  residuals <- x$residuals
  xName <- x$xName
  yName <- names(x$model)[1]
  distributionParams <- x$distributionParams
  coefficients <- x$coefficients


  hypDens <- function(residuals)
    dhyperb(residuals, param = distributionParams)
  logHypDens <- function(residuals)
    log(dhyperb(residuals, param = distributionParams))

  par(mfrow = c(2, 2), ...)

  plot(fitted.values, residuals, xlab = "Fitted Values",
       ylab = "Residuals", main = plotTitles[1])
  abline(h = 0, lty = 3)


  histData <- hist(residuals, breaks = breaks, right = FALSE,
                   plot = FALSE)
  breaks = histData$breaks
  empDens <- ifelse(!is.finite(log(histData$density)),
                    NA, histData$density)
  ymax <- 1.06 * max(hypDens(seq(min(breaks), max(breaks), 0.001)),
                     empDens, na.rm = TRUE)
    ##sde <- density(residuals, bw = 0.1)$y
    ##sden <- dhyperb(residuals, param = distributionParams)
    ##hist(residuals, freq = FALSE, breaks = breaks,
         ##ylim = c(0, max(sden, sde)),
         ##main = plotTitles[2], ...)
  hist(residuals, right = FALSE, freq = FALSE, ylim = c(0, ymax),
       main = plotTitles[2], breaks = breaks, ...)
  curve(dhyperb(x, param = distributionParams), add = TRUE,
        ylab = NULL, col = 2)


  logHist(residuals, breaks = breaks, include.lowest = TRUE,
          right = FALSE, main = plotTitles[3], ...)
  curve(log(dhyperb(x, param = distributionParams)),
        add = TRUE, ylab = NULL, xlab = NULL, col = 2)


  alphas <- (1:length(residuals) - 0.5)/length(residuals)
  quantiles <- qhyperb(alphas, param = distributionParams)
  qqplot(quantiles, residuals,
         main = plotTitles[4],
         xlab = "Hyperbolic quantiles",
         ylab = "Sample quantiles of residuals", ...)
  abline(0, 1)


}
