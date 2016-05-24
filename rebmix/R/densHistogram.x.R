.densHistogram.x <- function(k, x, x0, hx, cx, px)
{
  output <- .C("RdensHistogramX",
    k = as.integer(k),
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    x0 = as.double(x0),
    hx = as.double(hx),
    cx = as.integer(cx == .rebmix$Variables[2]),
    px = as.character(px),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densHistogram.x!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.x
