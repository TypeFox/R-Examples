writeJPEG <- function(image, target = raw(), quality = 0.7, bg = "white", color.space) {
  if (missing(color.space)) color.space <- attr(image, "color.space")
  if (inherits(target, "connection")) {
    r <- .Call("write_jpeg", image, raw(), quality, bg, color.space, PACKAGE="jpeg")
    writeBin(r, target)
    invisible(NULL)
  } else invisible(.Call("write_jpeg", image, if (is.raw(target)) target else path.expand(target), quality, bg, color.space, PACKAGE="jpeg"))
}
