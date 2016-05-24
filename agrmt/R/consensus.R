consensus <- function(V) {
# calculate consensus, according to Tatsle & Wierman (2007)
# arguments:   V = frequency vector

  # input validation: is this a frequency vector?
  if (min(V) < 0) stop("Error: negative values found in frequency vector.")
  if (all(round(V,0) == V) == FALSE) stop("Error: not all values are whole numbers; not a frequency vector; multiply by 100?")
 
  P <- V/sum(V) # normalizing the vector
  pos <- 1:length(V)
  mx <- mean(expand(V))
  dx <- max(pos) - min(pos)
  cns <- 1 + sum(P * log2(1-abs(pos-mx)/dx))
  # mathematically, consensus is not defined when all observations are the same
  # but by definition this is one:
  if (max(V)==sum(V)) cns <- 1
  return(cns)
}
