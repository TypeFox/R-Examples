## ipanel.smooth.R

ipanel.smooth <-
  function (x, y = NULL, pixs = 1, zmax = NULL,
            ztransf = function(x){x}, colramp = IDPcolorRamp,
            col = "black", lwd = 2, span = 2/3,
            iter = 3, ...)
  ## Author:  Rene Locher
  ## Version: 2007-02-07
{
  Image(x, y, pixs = pixs, ztransf = ztransf, colramp = colramp )
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = col, lwd = lwd, ...)
} ## ipanel.smooth
