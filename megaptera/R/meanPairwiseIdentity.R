meanPairwiseIdentity <- function(x){
  
  pairwiseIdentity <- function(x){
    sapply(seq_along(x), function(x, i) length(which(x[i] == x[-i])), x = x)
  }
  m <- apply(x, 2, pairwiseIdentity)
  m <- m/(nrow(x) - 1)
  colMeans(m)
}

myDist <- function(d){
  if (!inherits(d, "DNAbin")) stop("object not of class \"DNAbin\"")
  n <- ncol(d)
  d <- dist.dna(d, model = "N", pairwise.deletion = TRUE, 
                as.matrix = TRUE)
  d <- d/n
  d[d == 1] <- 0
  diag(d) <- NA
  return(d)
}


maxDistMPI <- function(a){
  
  if ( nrow(a) < 3 ){
    d <- 0
  } else {
    b <- as.raw(c(136, 40, 72, 24))
    percentInformation <- function(x, b){
      length(which(x %in% b))
    }
    m <- apply(a, 2, percentInformation, b)
    m <- m/nrow(a)
    a <- a[, m > .5]
    d <- myDist(a)
    d <- max(d, na.rm = TRUE)
  }
  d
}
