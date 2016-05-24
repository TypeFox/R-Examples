as.number <- function(x) {
  warn <- getOption("warn"); options(warn = -1)

  y <- try(as.numeric(x))

  options(warn = warn)

  if (any(is.na(y))) {
    as.character(x)
  }
  else {
    as.character(format(y))
  }
} ## as.number
