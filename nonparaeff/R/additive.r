additive <- function(base = NULL, frontier = NULL, noutput = 1){
  require(lpSolve)

  if(is.null(frontier))
    frontier <- base

  if(!is.null(base) & !is.null(frontier)){
    base <- as.matrix(base)
    frontier <- as.matrix(frontier)
  }

  if(ncol(base) != ncol(frontier))
    stop("Number of columns in base matrix and frontier matrix should be the same!")

  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)
  nf <- nrow(frontier)

  front.Y <- t(frontier[, 1:s])
  front.X <- t(frontier[, (s+1):(s+m)])
  base.Y <- t(base[, 1:s])
  base.X <- t(base[, (s+1):(s+m)])


  front.mat <- rbind(front.X, front.Y, 1)
  slack.mat <- rbind(diag(m+s) *c(rep(1, m), rep(-1, s)), rep(0, m +
  s))

  f.con <- cbind(front.mat, slack.mat)

  f.obj <- c(rep(0, nf), rep(1, m+s))
  f.dir <-  rep("==", m + s + 1)

  re <- data.frame(matrix(0, nrow = n, ncol = 1 + nf + m + s))
  names(re) <- c("eff",
                 paste("lambda", 1:nf, sep = ""),
                 paste("slack.x", 1:m, sep = ""),
                 paste("slack.y", 1:s, sep = ""))

  for(i in 1:n){
    f.rhs <- c(base.X[,i], front.Y[,i], 1)
    re.tmp <- lp("max", f.obj, f.con, f.dir, f.rhs)

    re[i,1] <- re.tmp$objval
    re[i, 2:(1+nf+m+s)] <- re.tmp$solution
  }
    
  return(re)
}
