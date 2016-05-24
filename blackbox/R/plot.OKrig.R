plot.OKrig <- function(x, digits = 4, which = 1, ...) {
  ## seul which=1 devrait marcher pour l'instant
  y <- x$y ## should work whatever the class of fitobject
  fitted.values <- x$fitted.values
  weights <- rep(1, length(fitted.values))
  #std.residuals <- (out$residuals * sqrt(out$weights))/out$shat.GCV
  if (any(which == 1)) {
    ##temp <- summary(out)
    #covariance <- cor(fitted.values * sqrt(weights), (y) * sqrt(weights))^2
    plot(fitted.values, y, ylab = "Estimates", xlab = "Predicted values", bty = "n", ...)
    abline(0, 1)
    hold <- par("usr")
    #        text(hold[1], hold[4], paste(" R^2 = ", format(round(100 *
    #            covariance, 2)), "%", sep = ""), cex = 0.8, adj = c(0, 1))
  }
}
