## rank correlation betwen distance on gradient (`g`)  and distance
## in response space (`y`)
## Similar to rankindex() in vegan but for analogue
rankDC <- function(g, y,
                   dc = c("chord", "bray", "euclidean", "chi.square", "gower"),
                   method = "spearman") {
  g <- as.data.frame(g)
  if (any(sapply(g, is.factor))) {
    message("'g' included factors: used Gower's distance")
    span <- distance(g, method = "mixed")
  } else {
    span <- distance(g, method = "euclidean")
  }
  span <- as.dist(span)
  y <- as.matrix(y)
  res <- numeric(length(dc))
  names(res) <- dc
  for (i in dc) {
    Y <- as.dist(distance(y, method = i))
    res[i] <- cor(span, Y, method = method)
  }
  class(res) <- "rankDC"
  res
}

print.rankDC <- function(x, ...) {
  attr(x, "class") <- NULL
  print.default(x)
}
