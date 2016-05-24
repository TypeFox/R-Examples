EF <- function (o, p) {
  N(o, p)
  EF <- 1 - sum((p - o)^2) / sum((o - mean(o))^2)
  EF
}

