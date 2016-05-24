baselineabif <- function(rfu, maxrfu = 1000){
  #
  # Check argument:
  #
  if(!is.numeric(rfu)) stop("numerical vector expected for rfu")
  #
  # Do not consider data above threshold maxrfu:
  #
  rfu[rfu >= maxrfu] <- NA
  #
  # Compute a kernel density estimate of data:
  #
  dst <- density(rfu, na.rm = TRUE)
  #
  # Choose as baseline the most common value:
  #
  baseline <- dst$x[which.max(dst$y)]
  return(baseline)
}
