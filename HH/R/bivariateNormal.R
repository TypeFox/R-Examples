bivariateNormal <-
  function(rho=0, layout=c(3,3), lwd=.2,
           angle=c(22.5, 67.5, 112.5, 337.5, 157.5, 292.5, 247.5, 202.5),
           col.regions=trellis.par.get("regions")$col, ...)
{
  x <- seq(-2, 2, length=33)
  y <- x
  fxy <- 1/(2*pi*sqrt(1-rho^2)) *
    exp(-.5/(1-rho^2) * outer(x^2, y^2, "+") - 2*rho*outer(x,y,"*"))
  n.angle <- length(angle)
  Angle <- rep(angle, rep(length(fxy), n.angle))
  Viewing.Angle <- ordered(Angle, angle)
  wireframe(rep(fxy, n.angle) ~ rep(x[row(fxy)], n.angle) * rep(y[col(fxy)], n.angle) |
            Viewing.Angle, r=rho, layout=layout, lwd=lwd, ...,
            panel = function(x, y, subscripts, z, angle, ...)
            {
              w <- unique(angle[subscripts])
              panel.wireframe(x=x, y=y, subscripts=subscripts, z=z,
                              screen = list(z = w, x = -60, y = 0), ...)
            },
            angle = Angle, ## this is how to pass down external element
            strip = strip.custom(strip.names = TRUE, style = "1", sep=": "),
            skip = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
            drape = TRUE, distance = 0.3,
            main = as.expression(substitute("Bivariate Normal, " * rho == r,
                                             c(alist(rho=rho), list(r=rho)))),
            xlab = list("x", cex = 0.6),
            ylab = list("y", cex = 0.6),
            zlab = list("f(x,y)", cex = 0.6))
}
