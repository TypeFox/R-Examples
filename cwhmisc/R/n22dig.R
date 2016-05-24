n22dig <- function(x, symm = TRUE) { ##
  n0 <- function(y) { ifelse(y == "100", " I", y) }
  s0 <- function(x) { if (nchar(x) >= 3 || substring(x,1,1)==" ") x <- substring(x, 2) else x }
  y <- format(round(100 * x))
  y1 <- n0(y)
  y2 <- matrix(sapply(y1,s0),nrow(x),ncol(x))
  if (is.matrix(x)) {
    if (symm) y2[row(x) < col(x)] <- "  "
    y2[y2=="I" ] <- " I"  # cosmetics
    y <- matrix(y2,nrow(x),ncol(x))  # should be unnecessary
  }  else {
    y2 <- n0(x)
}
  attr(y2,"names") <- NULL  ## names are a nuisance
  return(y2)
}
