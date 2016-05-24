pcplot <-
function(x, batch, y, alpha=0.35, ...) {
  ##require("MASS")
  p <- ncol(x)
  n <- nrow(x)
  xpr <- prcomp(x, scale. = FALSE)
  xp <- predict(xpr)[,1:2]
  plotc <- plotcomp(xp, groups=batch, y=y, alpha=alpha, ...)
  xm <- tapply(xp[,1], batch, mean)
  ym <- tapply(xp[,2], batch, mean)
  points(xm, ym, col=plotc$col, cex = 2, pch = 18)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}
