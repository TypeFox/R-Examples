readJPEG <- function(source, native=FALSE)
  .Call("read_jpeg", if (is.raw(source)) source else path.expand(source), native, PACKAGE="jpeg")
