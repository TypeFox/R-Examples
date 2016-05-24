direc.dea <-
  function(base = NULL, frontier = NULL, ngood = 1, nbad = 1)
{
  b.mat <- as.matrix(base)
  if(is.null(frontier))
    f.mat <- b.mat
  else
    f.mat <- as.matrix(frontier)

  if(ncol(b.mat) != ncol(f.mat))
    stop("Number of columns of the base matrix and the frontier matrix is not identical!")
  
  n <- nrow(b.mat)
  n1 <- nrow(f.mat)
  s <- ngood
  m <- nbad
  l <- ncol(b.mat) - s - m

  f.sign <- c(rep(-1, s), rep(1, m), rep(0, l))

  f.obj <- c(1, -1, rep(0, n1))
##  f.obj <- c(1, rep(0, n1))
  f.con1 <- t(f.mat)
  f.dir <- c(rep(">=", s), rep("==", m), rep("<=", l))

  re <- vector()
  
  for(i in 1:n){
    f.mat.n <- b.mat[i,]
    f.rhs <- f.mat.n
    f.con2 <- as.matrix(data.frame(f.sign*f.mat.n, -f.sign*f.mat.n))
##    f.con2 <- as.matrix(data.frame(f.sign*f.mat.n))
    f.con <- cbind(f.con2, f.con1)

    sol <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution
    re[i] <- sol[1] - sol[2]
##    re[i] <- sol[1]
  }
  return(re)
}
    
