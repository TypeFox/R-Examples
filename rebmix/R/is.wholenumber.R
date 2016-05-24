is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  warn <- getOption("warn"); options(warn = -1)

  x <- try(as.numeric(x))

  options(warn = warn)

  if (any(is.na(x))) {
    FALSE
  }
  else {
    is.numeric(x) && all(abs(x - round(x)) < tol)
  }
} ## is.wholenumber
