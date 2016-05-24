ci.plot <- function(x, y, se, level=0.95, ylim=NULL,
                    newplot=TRUE,
                    fun=gaussian()$linkinv,
                    dfun = gaussian()$mu.eta, ...) {
  z <- qnorm(1 - (1 - level) / 2)
  se <- dfun(y) * se
  y <- fun(y)
  u <- y + z * se
  l <- y - z * se
  if (is.null(ylim)) ylim <- range(c(u,l), na.rm=TRUE)
  if (newplot) plot(x, y, type="s", ylim=ylim, ...)
  else lines(x, y, type="s", ...)
  lines(x, u, type="s", lty=3, ...)
  lines(x, l, type="s", lty=3, ...)
}
