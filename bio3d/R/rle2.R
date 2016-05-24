rle2 <- function (x) {
  ## This is a modified version of base function "rle()"
  if (!is.vector(x) && !is.list(x)) 
    stop("'x' must be an atomic vector")
  n <- length(x)
  if (n == 0L) 
    return(structure(list(lengths = integer(), values = x, inds=integer), class = "rle2"))
  y <- x[-1L] != x[-n]
  i <- c(which(y | is.na(y)), n)
  structure(list(lengths = diff(c(0L, i)), values = x[i], inds=i), class = "rle2")
}

print.rle2 <- function(x, digits = getOption("digits"), prefix = "", ...) {
  if (is.null(digits)) 
    digits <- getOption("digits")
  cat("", "Run Length Encoding\n", "  lengths:", sep = prefix)
  utils::str(x$lengths)
  cat("", "  values :", sep = prefix)
  utils::str(x$values, digits.d = digits)
  cat("", "  indices:", sep = prefix)
  utils::str(x$inds, digits.d = digits)
  invisible(x)
}
