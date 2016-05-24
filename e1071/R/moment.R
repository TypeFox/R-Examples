moment <- function(x, order = 1, center = FALSE, absolute = FALSE,
		   na.rm = FALSE) {
  if (na.rm) 
    x <- x[!is.na(x)]
  if (center)
    x <- x - mean(x)
  if (absolute)
    x <- abs(x)
  sum(x ^ order) / length(x)
}
