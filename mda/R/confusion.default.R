confusion.default <-
  function (object, true, ...) 
{
  if (inherits(object, "data.frame")) 
    confusion.list(object, true)
  else {
    jt <- table(predicted=object, true)
    jd <- dimnames(jt)
    jn <- unlist(jd)
    ju <- jn[duplicated(jn)]
    j1 <- jd[[1]][!match(jd[[1]], ju, 0)]
    j2 <- jd[[2]][!match(jd[[2]], ju, 0)]
    jt <- jt[c(ju, j1), c(ju, j2), drop = FALSE]
    realjt <- unclass(jt[ju, ju, drop = FALSE])
    ntot <- sum(jt)
    mismatch <- (ntot - sum(realjt))/ntot
    structure(jt, error = (1 - sum(diag(realjt))/sum(realjt)), 
              mismatch = if (mismatch > 0) 
              mismatch
              else NULL)
  }
}

