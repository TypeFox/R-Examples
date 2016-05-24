print.audioInstance <- function(x, ...) {
  kind <- c("player", "recorder")[.Call("audio_instance_type", x, PACKAGE="audio")]
  info <- paste(" Audio ", kind, " instance ", sprintf('%x',.Call("audio_instance_address", x, PACKAGE="audio")), " of ", .Call("audio_driver_descr", x, PACKAGE="audio"), " (", .Call("audio_driver_name", x, PACKAGE="audio"), ").\n", sep = '')
  cat(info)
  invisible(info)
}

print.audioSample <- function(x, ...) {
  kind <- if (is.null(dim(x)) || dim(x)[1] != 2) 'mono' else 'stereo'
  bits <- attr(x, "bits", TRUE)
  bits <- if (is.null(bits)) '' else paste(", ", bits, "-bits", sep='')
  cat("sample rate: ", attr(x,"rate"), "Hz, ", kind, bits, "\n", sep='')
  attributes(x) <- NULL
  print(x)
}
