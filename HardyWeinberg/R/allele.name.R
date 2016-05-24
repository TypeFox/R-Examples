allele.name <- function(x,nr=1) {
  # extracts the name of the first (or the second) allele from the names of a vector of genotype counts.
  lx <- length(x)
  vec <- NULL
  for (i in 1:lx) {
    fa <- substr(x[i],nr,nr)
    vec <- c(vec,fa)
  }
  return(vec)
}
