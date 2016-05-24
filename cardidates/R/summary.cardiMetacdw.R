summary.cardiMetacdw <- function(object, file="", ...) {
  if (file=="") {
    print(object$metares)
  } else {
    write.table(object$metares, file=file, ...)
  }
}