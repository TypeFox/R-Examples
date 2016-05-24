dea <- function(base = NULL, frontier = NULL,
                noutput = 1, orientation=1, rts = 1, onlytheta = FALSE) {

  require(lpSolve)
  require(gdata)
  require(Hmisc)
  require(rms)

  if(is.null(frontier))
    frontier <- base

  if(!is.null(base) & !is.null(frontier)){
    if(is.null(nrow(base)))
      base <- matrix(base, nrow = 1)
    else
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
  if(n == 1){
    base.Y <- matrix(base[, 1:s], ncol = 1)
    base.X <- matrix(base[, (s+1):(s+m)], ncol = 1)
  }
  else{
    base.Y <- t(base[, 1:s])
    base.X <- t(base[, (s+1):(s+m)])
  }

  if(rts == 1){
    front.mat <- rbind(front.X, front.Y)
    slack.mat <- diag(m+s) *c(rep(1, m), rep(-1, s))
  }
  else{
    front.mat <- rbind(front.X, front.Y, 1)
    slack.mat <- rbind(diag(m+s) *c(rep(1, m), rep(-1, s)), rep(0, m + s))
  }
  
  first.front.mat <- front.mat
  second.front.mat <- cbind(front.mat, slack.mat)

  first.obj <- c(1, rep(0, nf))
  second.obj <- c(rep(0, nf), rep(1, m+s))

  if(onlytheta == FALSE){
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + nf + m + s))
    names(re) <- c("eff",
                   paste("lambda", 1:nf, sep = ""),
                   paste("slack.x", 1:m, sep = ""),
                   paste("slack.y", 1:s, sep = ""))
  }
  else{
    re <- data.frame(matrix(0, nrow = n, ncol = 1))
    names(re) <- c("eff")
  }
  
  for(i in 1:n){
    f.dir <- c(rep("<=", m), rep(">=", s))
    if(orientation == 1){
      f.con1 <- c(-base.X[, i], rep(0, s))
      f.rhs <- c(rep(0, m), base.Y[, i])
    }
    else{
      f.con1 <- c(rep(0, m), -base.Y[, i])
      f.rhs <- c(base.X[,i], rep(0, s))
    }

    if(rts == 2){
      f.con1 <- c(f.con1, 0)
      f.rhs <- c(f.rhs, 1)
      f.dir <- c(f.dir, "==")
    }

    f.con <- cbind(f.con1, first.front.mat)

    ## first step
    if(orientation == 1)
      tmp.re <- lp("min", first.obj, f.con, f.dir, f.rhs)
    else
      tmp.re <- lp("max", first.obj, f.con, f.dir, f.rhs)

    if(onlytheta == FALSE){
      re[i,1:(1+nf)] <- tmp.re$solution[1:(1+nf)]
    }
    else{
      re[i, 1] <- tmp.re$solution[1]
    }

    if(onlytheta == FALSE){
      ## second step
      f.obj <- c(rep(0, nf), rep(1, m+s))
      if(rts == 1){
        f.dir <- rep("==", m+s)
        f.rhs <- c(base.X[,i], base.Y[,i])
        if(orientation == 1)
          f.rhs <- f.rhs * c(rep(re[i,1], m), rep(1, s))
        else
          f.rhs <- f.rhs * c(rep(1, m), rep(re[i,1], s))
      }
      else{
        f.dir <- rep("==", m+s+1)
        f.rhs <- c(base.X[,i], base.Y[,i], 1)
        if(orientation == 1)
          f.rhs <- f.rhs * c(rep(re[i,1], m), rep(1, s), 1)
        else
          f.rhs <- f.rhs * c(rep(1, m), rep(re[i,1], s), 1)
      }

      f.con <- second.front.mat
      
      tmp.re <- lp("max", second.obj, f.con, f.dir, f.rhs)$solution
      re[i, (2+nf):(1+nf+m+s)] <- tmp.re[(nf+1):(nf+m+s)]
    }

  }
  return(re)
}
