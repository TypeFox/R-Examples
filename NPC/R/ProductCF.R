ProductCF <- function (p, B) {
  ## Fisher's combining function
  -2*sum(log(open01(p, B)))
}
