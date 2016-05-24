cont.table <- function (k1, n1, k2, n2, as.list=NA) {
  l <- max(length(k1), length(n1), length(k2), length(n2))
  if (missing(as.list)) as.list <- if (l > 1) TRUE else FALSE
  if (l > 1 && !as.list) stop("k1, n1, k2, n2 must be single numbers (as.list=FALSE)")

  # ensure that all input vectors have the correct length
  if (length(k1) != l) k1 <- rep(k1, length.out=l)
  if (length(n1) != l) n1 <- rep(n1, length.out=l)
  if (length(k2) != l) k2 <- rep(k2, length.out=l)
  if (length(n2) != l) n2 <- rep(n2, length.out=l)
  
  # sanity checks
  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  if (any(k1 + k2 <= 0)) stop("either k1 or k2 must be non-zero")

  # construct list of 2x2 contingency tables
  table.list <- lapply(1:l, function (i) matrix(c(k1[i], n1[i]-k1[i], k2[i], n2[i]-k2[i]), nrow=2, byrow=FALSE))
  
  if (as.list) table.list else table.list[[1]]
}
