play <- function(x, ...) UseMethod("play")
pause <- function(x, ...) UseMethod("pause")
resume <- function(x, ...) UseMethod("resume")
rewind <- function(x, ...) UseMethod("rewind")
wait <- function(x, ...) UseMethod("wait")

record <- function(where, rate, channels) {
  if (missing(rate)) {
    rate <- attr(where, "rate", TRUE)
    if (is.null(rate)) rate <- 44100
  }
  if (missing(channels))
    channels <- if (is.null(dim(where))) 2 else dim(where)[1]
  channels <- as.integer(channels)
  if (length(channels) != 1 || (channels != 1 && channels != 2))
    stop("channels must be 1 (mono) or 2 (stereo)")
  if (length(where) == 1) where <- if (channels == 2) matrix(NA_real_, 2, where) else rep(NA_real_, where)
  a <- .Call("audio_recorder", where, as.double(rate), as.integer(channels), PACKAGE="audio")
  .Call("audio_start", a, PACKAGE="audio")
  invisible(a)
}

pause.audioInstance <- function(x, ...)
  invisible(.Call("audio_pause", x, PACKAGE="audio"))

resume.audioInstance <- function(x, ...)
  invisible(.Call("audio_resume", x, PACKAGE="audio"))

rewind.audioInstance <- function(x, ...)
  invisible(.Call("audio_rewind", x, PACKAGE="audio"))

close.audioInstance <- function(con, ...)
  invisible(.Call("audio_close", con, PACKAGE="audio"))

wait.audioInstance <- function(x, timeout=NA, ...)
  invisible(.Call("audio_wait", x, if (any(is.na(timeout))) -1 else as.double(timeout), PACKAGE="audio"))

wait.default <- function(x, timeout, ...) {
  if (missing(timeout))
    timeout <- if (is.numeric(x)) x else NA
  invisible(.Call("audio_wait", NULL, if(any(is.na(timeout))) -1 else as.double(timeout), PACKAGE="audio"))
}

play.default <- function(x, rate=44100, ...) {
  a <- .Call("audio_player", x, rate, PACKAGE="audio")
  .Call("audio_start", a, PACKAGE="audio")
  invisible(a)
}

play.Sample <- function(x, ...) play(x$sound, x$rate)

play.audioSample <- function(x, rate, ...) {
  if (missing(rate)) rate <- attr(x, "rate", TRUE)
  play.default(x, rate, ...)
}

`[.audioSample` <- function(x, ..., drop = FALSE) {
  y <- NextMethod("[")
  attr(y, "rate") <- attr(x, "rate", TRUE)
  attr(y, "bits") <- attr(x, "bits", TRUE)
  class(y) <- class(x)
  y
}

play.audioInstance <- function(x, ...) stop("you cannot play an audio instance - try play(a$data) if a is a recorded instance")

`$.audioInstance` <- function(x, name) if (isTRUE(name == "data")) .Call("audio_instance_source", x, PACKAGE="audio") else NULL

`$.audioSample` <- function(x, name) attr(x, name)
`$<-.audioSample` <- function(x, name, value) .Primitive("attr<-")

audioSample <- function(x, rate=44100, bits=16, clip = TRUE) {
  if (!is.null(dim(x)) && dim(x)[1] != 1 && dim(x)[1] != 2)
    stop("invalid dimensions, audio samples must be either vectors or matrices with one (mono) or two (stereo) rows")
  if (is.integer(x)) {
    if (isTRUE(bits == 16)) x <- x / 32767.0 else if (isTRUE(bits == 8)) x <- x / 127.0 else stop("invalid sample size, must be 8 or 16 bits")
  }
  if (clip) {
    x[x > 1] <- 1
    x[x < -1] <- -1
  }
  attr(x, "rate") <- rate
  attr(x, "bits") <- as.integer(bits)
  class(x) <- "audioSample"
  x
}

as.audioSample <- function(x, ...) UseMethod("as.audioSample")

as.audioSample.default <- function(x, rate, bits, clip, ...) {
  if (missing(rate)) rate <- 44100
  if (missing(bits)) bits <- 16L
  if (missing(clip)) clip <- TRUE
  audioSample(x, rate, bits, clip)
}

# "sound" compatibility functions

as.audioSample.Sample <- function(x, ...)
  audioSample(x$sample, x$rate, x$bits)

as.Sample.audioSample <- function(x, ...) {
  y <- x
  attributes(y) <- NULL
  # we are not using sound's constructor, because we don't want to depend on it,
  # but it may prove to be dangerous if sound ever changes the format
  structure(list(sample = y, rate = x$rate, bits = x$bits), class="Sample")
}
