.estimateCov <- function(y1, y2) {
  y1.num <- nrow(y1)
  y2.num <- nrow(y2)
  gene.num <- ncol(y1)

  y1.mean        <- apply(y1, 2, mean)
  names(y1.mean) <- NULL
  y1.scal        <- y1 - matrix(rep(y1.mean, y1.num*gene.num), y1.num, gene.num, byrow=T)
  s1             <- cov(y1.scal)

  y2.mean        <- apply(y2, 2, mean)
  names(y2.mean) <- NULL
  y2.scal        <- y2 - matrix(rep(y2.mean, y2.num*gene.num), y2.num, gene.num, byrow=T)
  s2             <- cov(y2.scal)

  cov <- list(s1=s1, s2=s2)
  cov$s <- .estimateCovPool(y1.num, y2.num, s1, s2)
  return(cov)
}

.estimateCovPool <- function(y1.num, y2.num, s1, s2) {
  (s1*(y1.num-1)+ s2*(y2.num-1)) / (y1.num+y2.num-2)
}
