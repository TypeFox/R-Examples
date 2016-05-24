NormalCF <- function (p, B) {
  ## Liptak combining function
  -sum(stats::qnorm(open01(p, B)))
}
