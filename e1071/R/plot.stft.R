plot.stft <- function (x, col = gray (63:0/63), ...)
  {
    x <- x$values
    image(x=1:dim(x)[1], y=1:dim(x)[2], z=x, col=col, ...)
}
