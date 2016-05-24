mat.ip.l2d.kga <-
function(x)  {
  lot = x[, ncol(x)]
  x = x[, -ncol(x)]
  p = ncol(x)
  varL = by(x, lot, var)
  nb.groups = nlevels(lot)
  
  # Computation of the smoothing bandwidth matrices
  nbL<-by(x,INDICES=lot,FUN=nrow);
  wL<-bandwidth.parameter(p,nbL)
  # Multiplication of the variance by the window parameter
  varLwL<-varL
  for (i in 1:nb.groups)
    {varLwL[[i]]<-varL[[i]]*(wL[[i]]^2)}
  
  # Computation of the inner-products
  W = diag(0, nrow = nb.groups)
  dimnames(W) = list(levels(lot), levels(lot))
  xi = x[lot == levels(lot)[1], ]
  W[1, 1] = l2d.kgw(xi, varLwL[[1]], xi, varLwL[[1]])
  for (i in 2:nb.groups)  
    {xi = x[lot == levels(lot)[i], ]
    W[i, i] = l2d.kgw(xi, varLwL[[i]], xi, varLwL[[i]])
    for (j in 1:(i-1))  
      {xj = x[lot == levels(lot)[j], ]
      W[i, j] = W[j, i] = l2d.kgw(xi,varLwL[[i]],xj,varLwL[[j]])
    }
  }
  W
}
