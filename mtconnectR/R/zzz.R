.onLoad <- function(libname, pkgname) { # nocov start
  options(stringsAsFactors = FALSE)
  options(digits.secs = 6)
  options(scipen = 4)
  options(digits = 10)
  data.table::setNumericRounding(0L)
  Sys.setlocale('LC_ALL','C')
} # nocov end

Sys.setenv("TZ" = "UTC") # nocov
