ddf <- function(base = NULL, frontier = NULL,
                 noutput = 1, direction = NULL)
{
  # format of direction : (g_x's, g_y's)
  
  require(lpSolve)

  if(mode(direction) != "numeric" | is.null(direction))
    stop("direction must be a vector or NULL!!")
  
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

  if(is.null(direction))
    direction <- rep(1, m+s)

  gx <- direction[1:m]
  gy <- direction[(m+1):(m+s)]

  front.Y <- t(frontier[, 1:s])
  front.X <- t(frontier[, (s+1):(s+m)])
  base.Y <- t(base[, 1:s])
  base.X <- t(base[, (s+1):(s+m)])

  front.mat <- rbind(front.X, front.Y, 1)
  slack.mat <- rbind(diag(m+s) *c(rep(1, m), rep(-1, s)), rep(0, m + s))


  first.obj <- c(1, rep(0, nf))
  second.obj <- c(rep(0, nf), rep(1, m+s))

  re <- data.frame(matrix(0, nrow = n, ncol = 1 + nf + m + s))
  names(re) <- c("eff",
                 paste("lambda", 1:nf, sep = ""),
                 paste("slack.x", 1:m, sep = ""),
                 paste("slack.y", 1:s, sep = ""))

  for(i in 1:n){
    ## first step
    f.obj <- c(1, rep(0, nf))
    f.dir <- c(rep("<=", m), rep(">=", s), "==")

    f.con1 <- c(gx, -gy, 0)
    f.con <- cbind(f.con1, front.mat)
    f.rhs <- c(base.X[, i], base.Y[, i], 1)

    re.tmp <- lp("max", f.obj, f.con, f.dir, f.rhs)
    theta <- re.tmp$solution[1]
    lambda <- re.tmp$solution[2:(1+nf)]
    sx <- -front.X %*% as.matrix(lambda, ncol = 1) - theta*gx +
    base.X[,i]
    sy <- front.Y %*% as.matrix(lambda, ncol = 1) - theta*gy -
    base.Y[,i]

    sx[sx < 1e-10] <- 0
    sy[sy < 1e-10] <- 0


    re[i, 1] <- theta
    re[i, 2:(1+nf+m+s)] <- c(lambda, sx, sy)
  }

  return(re)
}
    
    

    
  
