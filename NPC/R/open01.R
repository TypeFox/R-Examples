open01 <- function (p, B) {
  ## map pvalue to open interval (0, 1)
  (p + 1/(2*B))/(1 + 1/B)
}
