# CUSTOM DISTANCE CALCULATOR
# package: megaptera
# called by: stepE
# author: Christoph Heibl
# last update: 2014-09-24


myDist <- function(d){
  if (!inherits(d, "DNAbin")) stop("object not of class \"DNAbin\"")
  n <- ncol(d)
  d <- dist.dna(d, model = "N", pairwise.deletion = TRUE, 
                as.matrix = TRUE)
  d <- d/n
  d[d == 1] <- 0
  return(d)
}