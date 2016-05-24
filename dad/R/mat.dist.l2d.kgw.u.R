mat.dist.l2d.kgw.u <-
function(x, varwL)  {
  lot = x[, ncol(x)]
  x = x[, -ncol(x)]
  distances = diag(0, nrow = nlevels(lot))
  dimnames(distances) = list(levels(lot), levels(lot))
  for (i in 2:nlevels(lot))  for (j in 1:(i-1))  {
    i.lot = which(lot == levels(lot)[i])
    j.lot = which(lot == levels(lot)[j])
    distances[i, j] = distances[j, i] = dist.l2d.kgw.u(x[i.lot], varwL[[i]], x[j.lot], varwL[[j]])
  }
  as.dist(distances)
}
