.densHistogram.xy <- function(k, x, y, x0, y0, hx, hy, cx, cy, px, py)
{
  output <- .C("RdensHistogramXY",
    k = as.integer(k),
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    x0 = as.double(x0),
    y0 = as.double(y0),
    hx = as.double(hx),
    hy = as.double(hy),
    cx = as.integer(cx == .rebmix$Variables[2]),
    cy = as.integer(cy == .rebmix$Variables[2]),
    px = as.character(px),
    py = as.character(py),        
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densHistogram.xy!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k
  length(output$z) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.xy
