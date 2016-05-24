#.First.lib <- function(lib, pkg)
#        library.dynam("cts", pkg, lib)
spectrum <- function(object, ...)
    UseMethod("spectrum")

"spectrum.car" <-
function (object, frmult=1, n.freq, plot.it = TRUE, na.action = na.fail, ...) 
{
  x <- object
  ## must be a result of an AR fit
#  cn <- match(c("frequency", "spectrum"), names(x))
  cn <- match(c("b","order","scale"),names(x))
  if (any(is.na(cn))) 
    stop("object must be a car() fit")
  b <- x$b
  order <- x$order
  scale <- x$scale
  if(missing(n.freq)) n.freq <- 500
  Z <- .Fortran("cspec",
           as.double(b),
           as.integer(order),
           as.double(scale),
           as.double(frmult),
           as.integer(n.freq),
           freq=double(n.freq),
           spec=double(n.freq),
           package="cts")
  spg.out <- list(freq = Z$freq, spec = Z$spec, method = paste("CAR (",order, ") spectrum ", sep = ""))
  class(spg.out) <- "spec.car"
  if (plot.it) {
    plotSpecCar(spg.out, ci = 0, ...)
    return(invisible(spg.out))
  }
  else return(spg.out)
}
