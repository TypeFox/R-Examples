block.smooth.periodogram = function(y, x = NULL, N = NULL, S = NULL, p = 0.25, spar.freq = 0, spar.time = 0, theta = 0, phi = 0, xlim = NULL, ylim = NULL, zlim = NULL, ylab = "Time", palette.col =  NULL){
  T. = length(y)
  if(is.null(N)){N = trunc(T.^0.8)}
  if(is.null(S)){S = trunc(p * N)}
  M  = trunc((T.-N)/S+1)
  aux = matrix(NA, ncol = M, nrow = trunc(N/2))
  for(j in 1:M){
    aux[,j] = periodogram(y[(S*(j-1)+1):(S*(j - 1)+N)], plot = FALSE)$periodogram
  }
  lambda = periodogram(y[(S*(j-1)+1):(S*(j - 1)+N)], plot = FALSE)$lambda
  
  aux2 = aux
  for(j in 1:M){
    aux2[,j] = smooth.spline(aux[,j], spar = spar.freq)$y
  }
  
  aux3 = aux
  for(i in 1:(dim(aux)[1])){
    aux3[i,] = smooth.spline(aux2[i,], spar = spar.time)$y
  }
  aux = aux3
  
  nrz = nrow(aux)
  ncz = ncol(aux)
  if(is.null(palette.col)){palette.col = c("green", "lightgreen", "yellow", "orange", "darkred")}
  jet.colors = colorRampPalette(palette.col)
  nbcol = 100
  color = jet.colors(nbcol)
  z = aux
  zfacet = z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol = cut(zfacet, nbcol)
  
  j = 1:M
  t = S * (j - 1) + N/2
  if(is.null(x)){t = t}
  else{t = x[t]}
  if(is.null(xlim)){xlim = range(lambda)}
  if(is.null(ylim)){ylim = range(t)}
  if(is.null(zlim)){zlim = range(aux, na.rm = TRUE)}
  persp(x = lambda, y = t, z = aux, theta = theta, phi = phi, col = color[facetcol], zlab = "Smooth Periodogram", xlab = "Frequency", ylab = ylab, expand = 0.5, ticktype = "detailed", ylim = ylim, zlim = zlim, xlim = xlim)
}