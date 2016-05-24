centers.interval <-
function(sym.data) {
  idn <- all(sym.data$sym.var.types == sym.data$sym.var.types[1])
  if (idn == FALSE) 
    stop("All variables have to be of the same type")
  
  if ((sym.data$sym.var.types[1] != "$I")) 
    stop("Variables have to be continuos or Interval")
  else
    nn <- sym.data$N
  mm <- sym.data$M
  centers <- matrix(0, nn, mm)
  centers <- as.data.frame(centers)
  rownames(centers) <- sym.data$sym.obj.names
  colnames(centers) <- sym.data$sym.var.names
  for (i in 1:nn) 
    for (j in 1:mm) 
      centers[i, j] <- (sym.var(sym.data,j)$var.data.vector[i,1] + sym.var(sym.data, j)$var.data.vector[i, 2])/2
  
  return(centers)
}
