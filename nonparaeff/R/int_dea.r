int.dea <- function(base = NULL, frontier = NULL,
                 noutput = 1, intinput = 1, orientation=1,
                 epsilon = 1e-06) {
  ## data.frame(ys, xi, xc)

  require(lpSolve)

  if(intinput < 1){
    stop("The number of integer input must be larger than 1.")
  }
  
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
  p <- intinput ## number of integer input
  q <- m - p    ## number of continuout input
  n <- nrow(base)
  nf <- nrow(frontier)

  front.Y <- t(frontier[, 1:s])
  front.Xi <- t(frontier[, (s+1):(s+p)])
  if(q != 0)
    front.Xc <- t(frontier[, (s+p+1):(s+p+q)])
  base.Y <- t(base[, 1:s])
  base.Xi <- t(base[, (s+1):(s+p)])
  if(q != 0)
    base.Xc <- t(base[, (s+p+1):(s+p+q)])

  if(q != 0)
    slack.mat <-
      rbind(cbind(-diag(s), matrix(0, s, q + p + p + p)),
            cbind(matrix(0, q, s), diag(q), matrix(0, q, p + p + p)),
            cbind(matrix(0, p, s + q), diag(p), matrix(0, p, p),
                  -diag(p)),
            cbind(matrix(0, p, s + q + p), diag(p), diag(p)))
  else
    slack.mat <-
      rbind(cbind(-diag(s), matrix(0, s, p + p + p)),
            cbind(matrix(0, p, s), diag(p), matrix(0, p, p),
                  -diag(p)),
            cbind(matrix(0, p, s + p), diag(p), diag(p)))

  if(q != 0)
    con1 <-
      rbind(front.Y, front.Xc, front.Xi, matrix(0, p, nf))
  else
    con1 <-
      rbind(front.Y, front.Xi, matrix(0, p, nf))
  
  f.con1 <- cbind(con1, slack.mat)

  if(q != 0){
    f.obj <- c(1, rep(0, nf), -epsilon * rep(1, s + m + p), rep(0, p))
    f.dir <- rep("==", s + p + m)
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + nf + s + q + p + p
    + p ))
    colnames(re) <- c("theta", paste("lambda", 1:nf, sep = ""),
                      paste("s^+_", 1:s, sep = ""),
                      paste("s^-_", 1:m, sep = ""),
                      paste("s^I_", 1:p, sep = ""),
                      paste("til.x", 1:p, sep = ""))
  }
  else{
    f.obj <- c(1, rep(0, nf), -epsilon*rep(1, s + p + p), rep(0, p))
    f.dir <- rep("==", s + p + p)
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + nf + s + p + p +
    p))
    colnames(re) <- c("theta", paste("lambda", 1:nf, sep = ""),
                      paste("s^+_", 1:s, sep = ""),
                      paste("s^-_", 1:p, sep = ""),
                      paste("s^I_", 1:p, sep = ""),
                      paste("til.x", 1:p, sep = ""))
  }


  for(i in 1:n){
    if(q != 0)
      f.con2 <- c(rep(0, s), -base.Xc[, i], rep(0, p), -base.Xi[, i])
    else
      f.con2 <- c(rep(0, s), rep(0, p), -base.Xi[, i])
    
    f.con <- cbind(f.con2, f.con1)

    f.rhs <- c(base.Y[, i], rep(0, m + p))

    if(q != 0)
      sol <- lp("min", f.obj, f.con, f.dir, f.rhs,
                int.vec = (1+ nf + s + m + p + 1):(1+ nf + s + m + p +
                                                   p))
    else
      sol <- lp("min", f.obj, f.con, f.dir, f.rhs,
                int.vec = (1 + nf + s + p + p + 1):(1 + nf + s + p + p
                                                    + p))
    re[i,] <- sol$solution
  }

  return(re)
}
