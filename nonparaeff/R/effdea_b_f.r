effdea.b.f <- function(base = NULL, frontier = NULL, noutput = 1,
                       orientation=1, rts = 1, convhull = TRUE) {
  require(lpSolve)
  require(geometry)

  if(is.null(frontier))
    frontier <- base

  if(!is.null(base) & !is.null(frontier)){
    base <- as.matrix(base)
    frontier <- as.matrix(frontier)
  }

  if(ncol(base) != ncol(frontier))
    stop("Number of columns in base matrix and frontier matrix should be the same!")

  s <- noutput
  m <- dim(base)[2] - s

  n <- dim(frontier)[1]

  if(convhull == TRUE)
    idx <- sort(unique(as.vector(convhulln(frontier))))
  else
    idx <- 1:nrow(frontier)

  frontier <- frontier[idx,]

  nf <- nrow(frontier)

  base.Y <- t(base[, 1:s])
  base.X <- t(base[, (s+1):(s+m)])
  front.Y <- t(frontier[, 1:s])
  front.X <- t(frontier[, (s+1):(s+m)])

  re <- data.frame(matrix(0, n, 1))
  colnames(re) <- "eff"
  
  f.obj <- c(1, rep(0, nf))
  f.con1 <- rbind(front.X, front.Y)

  for(i in 1:n){
    if(rts == 1){
      f.dir <- c(rep("<=", m), rep(">=", s))
      if(orientation == 1){
        f.con2 <- c(-base.X[,i], rep(0, s))
        f.con <- cbind(f.con2, f.con1)
        f.rhs <- c(rep(0, m), base.Y[,i])
        re[i,1] <- lp("min", f.obj, f.con, f.dir, f.rhs)$objval
      }
      else{
        f.con2 <- c(rep(0, m), -base.Y[,i])
        f.con <- cbind(f.con2, f.con1)
        f.rhs <- c(base.X[,i], rep(0, s))
        re[i,1] <- lp("max", f.obj, f.con, f.dir, f.rhs)$objval
      }
    }
    if(rts == 2){
      f.dir <- c(rep("<=", m), rep(">=", s), "==")
      if(orientation == 1){
        f.con2 <- c(-base.X[,i], rep(0, s))
        f.con <- rbind(cbind(f.con2, f.con1), c(0, rep(1, nf)))
        f.rhs <- c(rep(0, m), base.Y[,i], 1)
        re[i,1] <- lp("min", f.obj, f.con, f.dir, f.rhs)$objval
      }
      else{
        f.con2 <- c(rep(0, m), -base.Y[,i])
        f.con <- rbind(cbind(f.con2, f.con1), c(0, rep(1, nf)))
        f.rhs <- c(base.X[,i], rep(0, s), 1)
        re[i,1] <- lp("max", f.obj, f.con, f.dir, f.rhs)$objval
      }
    }
  }
  return(re)
}
