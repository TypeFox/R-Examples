plot.fmres <- function (x, type = c("image", "fts", "contour", "filled.contour", "persp"), xlab = "Year", 
    ylab = "Age", zlab = "Residual", ...) 
{
    type <- match.arg(type)
    switch(type, 
      image = image(x$x, x$y, x$z, ylab = ylab, xlab = xlab, col = diverge_hcl(25, c=360, l=c(10,90)), ...), 
      fts = plot(fts(x$y, t(x$z), start = x$x[1], xname = xlab, yname = "Residuals"), ...), 
      contour = contour(x$x, x$y, x$z, ylab = ylab, xlab = xlab, ...), 
      filled.contour = filled.contour(x$x, x$y, x$z, ylab = ylab, xlab = xlab, ...), 
      persp = persp(x$x, x$y, x$z, ylab = ylab, xlab = xlab, zlab = zlab, ...)
    )
}
