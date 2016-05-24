chisq <- function (k1, n1, k2, n2, correct=TRUE, one.sided=FALSE) {
  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  if (any(k1 + k2 <= 0)) stop("either k1 or k2 must be non-zero")

  l <- max(length(k1), length(n1), length(k2), length(n2)) # ensure that all vectors have the same length
  if (length(k1) < l) k1 <- rep(k1, length.out=l)
  if (length(n1) < l) n1 <- rep(n1, length.out=l)
  if (length(k2) < l) k2 <- rep(k2, length.out=l)
  if (length(n2) < l) n2 <- rep(n2, length.out=l)

  k1 <- as.numeric(k1)                  # force integer -> float conversion to avoid overflow in multiplication below
  n1 <- as.numeric(n1)
  k2 <- as.numeric(k2)
  n2 <- as.numeric(n2)
  
  O11 <- k1                             # construct "observed" contingency table
  O21 <- n1 - k1
  O12 <- k2
  O22 <- n2 - k2

  R1 <- O11 + O12                       # compute row/column sums and sample size
  R2 <- O21 + O22
  C1 <- n1
  C2 <- n2
  N <- n1 + n2

  ## common form for homogeneity test with Yates' correction (Evert 2004, p.82)
  term <- abs(O11 * O22 - O12 * O21)
  if (correct) term <- pmax(term - N/2, 0)
  X2 <- (N * term^2) / (R1 * R2 * C1 * C2)

  # approximate one-sided chi-squared statistic as signed root of X2 (-> standard normal distribution)
  if (one.sided) {                      
    X2 <- sign(k1/n1 - k2/n2) * sqrt(X2)
  }

  X2
}
