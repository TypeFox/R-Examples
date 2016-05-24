GRI <- function (o, p) {
  N(o, p)
  expr <- sqrt(mean(((p - o) / (p + o))^2))
  GRI <- (1 + expr) / (1 - expr)
  GRI
}

