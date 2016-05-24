reQ <-
function(x, dist, ties="random") {
  if(length(x) != sum(dist)) {stop("Length of x must be equal to sum of dist.")}
  Q.scores <- rep(1:length(dist), times=dist)
  if(is.matrix(x) | is.data.frame(x)) {
    res <- t(apply(x, 1, function(j) Q.scores[rank(j, ties.method=ties)]))
    colnames(res) <- paste(colnames(x), ".reQ", sep="")
    return(res)
  }
  res <- Q.scores[rank(x, ties.method=ties)]
  res
}
