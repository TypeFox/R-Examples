is.number <- function(x) {
  warn <- getOption("warn"); options(warn = -1)

  y <- try(as.numeric(x))

  options(warn = warn)

  if (any(is.na(y))) {
    FALSE
  }
  else {
    TRUE
  }
} ## is.number
