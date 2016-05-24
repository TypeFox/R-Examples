
rvsample <- function(x, size=1, jointly=TRUE, reject.na=FALSE)
{
  # NAME
  #   rvsample - Draw Samples from Random Vectors
  #
  xs <- sims(as.rv(x))
  ns <- nrow(xs)
  if (is.null(size) || is.na(size)) size <- ns
  if (jointly) {
    if (reject.na) {
      f <- function (x) any(is.na(x))
      is.na.xs <- apply(xs, 1, f)
      if (all(is.na.xs)) {
        s <- sample(ns, size=size, replace=TRUE, prob=is.na.xs)
      } else {
        s <- sample(ns, size=size, replace=TRUE, prob=!is.na.xs) 
      }
    } else {
      s <- sample(ns, size=size, replace=TRUE)
    }
    s <- xs[s,]
  } else {
    s <- apply(xs, 2, function (s) {
      if (all(nas <- is.na(s))) return(s)
      sample(s[!nas], size=size, replace=TRUE)
    })
  }
  return(s)
}
