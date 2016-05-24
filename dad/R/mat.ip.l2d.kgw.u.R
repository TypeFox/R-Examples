mat.ip.l2d.kgw.u <-
function(x, varwL)  {
  lot = x[, ncol(x)]
  x = x[, -ncol(x)]
  W = diag(0, nrow = nlevels(lot))
  dimnames(W) = list(levels(lot), levels(lot))
  xi = x[lot == levels(lot)[1]]
  W[1, 1] = l2d.kgw.u(xi, varwL[[1]], xi, varwL[[1]])
  for (i in 2:nlevels(lot))
    {xi = x[lot == levels(lot)[i]]
    W[i, i] = l2d.kgw.u(xi, varwL[[i]], xi, varwL[[i]])
    for (j in 1:(i-1))
      {xj = x[lot == levels(lot)[j]]
      W[i, j] = W[j, i] = l2d.kgw.u(xi, varwL[[i]], xj, varwL[[j]])
    }
  }
  W
}
