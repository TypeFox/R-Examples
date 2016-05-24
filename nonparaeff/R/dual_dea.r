dual.dea <- function(base = NULL, frontier = NULL,
                     noutput = 1, orientation=1, rts = 1) {

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

  front.Y <- as.matrix(frontier[, 1:s])
  front.X <- as.matrix(frontier[, (s+1):(s+m)])
  base.Y <- as.matrix(base[, 1:s])
  base.X <- as.matrix(base[, (s+1):(s+m)])

  re <- data.frame(matrix(0, nrow = n, ncol = 1 + m + s))
  names(re) <- c("objval", paste("u", 1:s, sep = ""),
                 paste("v", 1:m, sep = ""))
  if(rts == 2){
    re <- data.frame(matrix(0, nrow = n, ncol = 2 + m + s))
    names(re) <- c("objval", paste("u", 1:s, sep = ""),
                   paste("v", 1:m, sep = ""), "z")
  }

  for(i in 1:n){
    if(orientation == 1){
      f.obj <- c(rep(base.Y[i,], s), rep(0, m))
      f.con1 <- c(rep(0, s), base.X[i,])
      f.con2 <- cbind(base.Y, -base.X)
      f.con <- rbind(f.con1, f.con2)
      f.dir <- c("==", rep("<=", nf))
      f.rhs <- c(1, rep(0, nf))
      if(rts == 2){
        f.obj <- c(f.obj, -1, 1)
        z.tmp <- c(0, rep(1, nf))
        z <- cbind(-z.tmp, z.tmp)
        f.con <- cbind(f.con, z)
      }
      re.tmp <- lp("max", f.obj, f.con, f.dir, f.rhs)
      re.obj <- re.tmp$objval
      re.uvz.tmp <- re.tmp$solution
      if(rts == 1)
        re.uvz <- re.uvz.tmp
      else
        re.uvz <- c(re.uvz.tmp[1:(length(re.uvz.tmp)-2)],
                    re.uvz.tmp[length(re.uvz.tmp)-1] -
                    re.uvz.tmp[length(re.uvz.tmp)])
      re[i, ] <- c(re.obj, re.uvz)
    }

    if(orientation == 2){
      f.obj <- c(rep(0, s), rep(base.X[i,], m))
      f.con1 <- c(base.Y[i,], rep(0, m))
      f.con2 <- cbind(base.Y, -base.X)
      f.con <- rbind(f.con1, f.con2)
      f.dir <- c("==", rep("<=", nf))
      f.rhs <- c(1, rep(0, nf))
      if(rts == 2){
        f.obj <- c(f.obj, -1, 1)
        z.tmp <- c(0, rep(1, nf))
        z <- cbind(z.tmp, -z.tmp)
        f.con <- cbind(f.con, z)
      }
      re.tmp <- lp("min", f.obj, f.con, f.dir, f.rhs)
      re.obj <- re.tmp$objval
      re.uvz.tmp <- re.tmp$solution
      if(rts == 1)
        re.uvz <- re.uvz.tmp
      else
        re.uvz <- c(re.uvz.tmp[1:(length(re.uvz.tmp)-2)],
                    re.uvz.tmp[length(re.uvz.tmp)-1] -
                    re.uvz.tmp[length(re.uvz.tmp)])
      re[i, ] <- c(re.obj, re.uvz)
    }
  }
  return(re)
}
