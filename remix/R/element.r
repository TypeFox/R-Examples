## elements <- function(x) {
##   xx <- paste(x, "~ 1")
##   f <- as.formula(xx)
##   if (grepl("^Surv\\(", x) | grepl("^I\\(", x) | grepl("factor\\(", x))
##     res <- x
##   else
##     res <- sapply(f[[2]], function(y) deparse(y, 500))
##   if (res[1] == ":") {
##     res <- paste(res[-1], collapse = ":")
##   }
##   res[res != "cbind"]
## }

elements <- function(x) {
  xx <- paste(x, "~ 1")
  f <- as.formula(xx)
  if (grepl("^cbind\\(", x))
    res <- sapply(f[[2]], function(y) deparse(y, 500))
  else
    res <- x
  res[res != "cbind"]
}
