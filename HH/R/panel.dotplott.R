## S-Plus only.
## This will never be called from R, so it is not inside an if.R() construction
## identical to panel.dotplot except for the transpose argument.

panel.dotplott <-
function(x, y,
         pch = dot.symbol$pch, col = dot.symbol$col,
         cex = dot.symbol$cex, font = dot.symbol$font, ..., transpose=FALSE)
{
  ok <- !is.na(x) & !is.na(y)
  dot.symbol <- trellis.par.get("dot.symbol")
  dot.line <- trellis.par.get("dot.line")
  if (transpose)
    abline(v = unique(x[ok]),
           lwd = dot.line$lwd, lty = dot.line$lty, col = dot.line$col)
  else
    abline(h = unique(y[ok]),
           lwd = dot.line$lwd, lty = dot.line$lty, col = dot.line$col)
  points(x, y, pch = pch, col = col, cex = cex, font = font, ...)
}
