multipheno.T2 <- function(z) {
  z <- as.matrix(z)
  ic <- solve(cor(z, method = "spearman"))
  zic <- z %*% ic
  t2 <- sapply(1:nrow(z), function(idx) return(sum(zic[idx, , drop = TRUE] * z[idx, , drop = TRUE])))
  pval <- pchisq(t2, df = ncol(z), lower.tail = FALSE)
  return(list(T2 = t2, pval = pval))
}
