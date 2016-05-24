entropy <- function(V) {
# calculates the Shannon entropy, following Tatsle & Wierman (2007)
  V <- V/sum(V)    # normalizing the vector
  a <- log2(V)     # log2 separately, so I can filter out the -Inf
  # this is how Tatsle & Wierman handle it!
  a[a==-Inf] <- 0  # if log(0), then 0
  e <- -1 * sum(V * a)
  return(e)
}