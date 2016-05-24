`aspectSMC` <-
function(r, targvar) {
  nvar <- targvar
  m <- dim(r)[1]
  ind <- which(1:m!=nvar)
  rn <- as.vector(r[nvar,ind])
  rr <- r[ind,ind]
  b <- solve(rr,rn)
  g <- diag(m)
  g[nvar,ind] <- -b
  g[ind,nvar] <- -b
  g[ind,ind] <- outer(b,b)
  list(f = sum(b*rn), g = -g)
}

