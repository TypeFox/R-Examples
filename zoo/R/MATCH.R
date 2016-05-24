MATCH <- function(x, table, nomatch = NA, ...)
  UseMethod("MATCH")
  
MATCH.default <- function(x, table, nomatch = NA, ...)
  match(x, table, nomatch = nomatch, ...)

MATCH.timeDate <- function(x, table, nomatch = NA, ...) {
  match(as.POSIXct(x), as.POSIXct(table), nomatch = nomatch, ...)
}

MATCH.times <- function(x, table, nomatch = NA, units = "sec", eps = 1e-10, ...) {
 match(trunc(x, units, eps), trunc(table, units, eps), nomatch = nomatch, ...)
}

