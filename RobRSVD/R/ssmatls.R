ssmatls <-
function(n){
  h = 1/(n+1)
  # Why not 1/(n-1) since n is the number of points ?
  Qmat = mat.or.vec(n, n-1)
  for (j in 2:(n-1)){
    Qmat[j-1,j] = 1/h
    Qmat[j,j] = -2/h
    Qmat[j+1,j] = 1/h
  }
  Qmat = Qmat[,2:(n-1)]
  
  Rmat = mat.or.vec(n-1, n-1)
  for (i in 1:(n-2)){
    Rmat[i,i] = (2/3)*h
    Rmat[i,i+1] = h/6
    Rmat[i+1,i] = h/6
  }
  Rmat[n-1, n-1] = (2/3)*h
  Rmat = Rmat[2:(n-1), 2:(n-1)]
  
  y = Qmat %*% solve(Rmat) %*% t(Qmat)
  return(list(y=y, Qmat=Qmat, Rmat=Rmat, h=h))
}
