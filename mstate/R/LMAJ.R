LMAJ <- function(msdata, s, from, method=c("aalen", "greenwood"))
{
  tmat <- attr(msdata, "trans")
  if (is.null(tmat)) stop("msdata object should have a \"trans\" attribute")
  K <- nrow(tmat)
  if (any(is.na(match(from, 1:K)))) stop("from should be subset of 1:K with K number of states")
  xss <- xsect(msdata, s)
  infrom <- xss$id[xss$state %in% from]
  msdatas <- cutLMms(msdata, LM=s)
  msdatasfrom <- msdatas[msdatas$id %in% infrom, ]
  c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=msdatasfrom)
  msf0 <- msfit(c0, trans=tmat)
  pt0 <- probtrans(msf0, predt=s, method=method)[from]
  if (length(from) == 1)
    return(pt0[[1]])
  else {
    xsss <- xss[xss$state %in% from, ]
    xsss$state <- factor(xsss$state, levels=from)
    tbl <- table(xsss$state)
    p <- tbl / sum(tbl)
    varp <- (diag(p) - p %*% t(p)) / sum(tbl)
    # Relying on fact that all items from list have exactly the same size and structure
    # and that the sum of p equals 1; sorry for the double for-loop
    res <- tmp1 <- tmp2 <- 0
    for (j in 1:length(from)) {
      ptj <- pt0[[j]]
      res <- res + p[j] * ptj[, 1 + (1:K)]
      tmp2 <- tmp2 + p[j] * (1 - p[j]) * (ptj[, K+1 + (1:K)])^2
      for (k in 1:length(from)) {
        ptk <- pt0[[k]]
        tmp1 <- tmp1 + varp[j,k] * ptj[, 1 + (1:K)] * ptk[, 1 + (1:K)]
      }
    }
    ses <- sqrt(tmp1 + tmp2)
    # take ptj as template and insert estimates and SE's
    pt <- ptj
    pt[, 1 + (1:K)] <- res
    pt[, K+1 + (1:K)] <- ses
    return(pt)
  }
}
