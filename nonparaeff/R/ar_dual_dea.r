ar.dual.dea <-
  function(base = NULL, frontier = NULL,
           noutput = 1, orientation=1, rts = 1, ar.l = NULL,
           ar.r = NULL, ar.dir = NULL, dual = FALSE)
{
  require(lpSolve)

  if(is.null(ar.l) | !is.matrix(ar.l))
    stop("You should give the left hand side of AR restrictions as a form of a matrix!")
  if(is.null(ar.r))
    stop("You should give the right hand side of AR restrictions as a form of a matrix or vector!")
  if(nrow(ar.l) != length(as.vector(ar.r)))
    stop("The number of rows in the left hand side of AR restriction matrix should be the same as the number of rows in the right hand side of AR restrictions")
  if(is.null(ar.dir))
    stop("You should give (in)equality for AR restrictions")
  if(length(ar.dir) != length(as.vector(ar.r)))
    stop("The length of ar.dir and the length of ar.r are different.")
  if(length(ar.dir) != nrow(ar.l))
    stop("The length of ar.dir and the number of columns in the left hand side of AR restriction are different.")
  
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

  if(ncol(ar.l) != s + m)
    stop("The number of columns in AR restriction matrix should be the total sum o f the number of inputs and outputs!!")
  
  front.Y <- as.matrix(frontier[, 1:s])
  front.X <- as.matrix(frontier[, (s+1):(s+m)])
  base.Y <- as.matrix(base[, 1:s])
  base.X <- as.matrix(base[, (s+1):(s+m)])

  if(dual == FALSE){
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + m + s))
    names(re) <- c("objval", paste("u", 1:s, sep = ""),
                   paste("v", 1:m, sep = ""))

    if(rts == 2){
      re <- data.frame(matrix(0, nrow = n, ncol = 2 + m + s))
      names(re) <- c("objval", paste("u", 1:s, sep = ""),
                     paste("v", 1:m, sep = ""), "z")
    }
  }
  else{
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + nf + nrow(ar.l) +
                            m + s))
    names(re) <- c("eff", paste("lambda", 1:nf, sep = ""),
                   paste("q", 1:nrow(ar.l), sep = ""),
                   paste("slack.x", 1:m, sep = ""),
                   paste("slack.y", 1:s, sep = ""))
  }
  

  for(i in 1:n){
    if(orientation == 1){
      f.obj <- c(rep(base.Y[i,], s), rep(0, m))
      f.con1 <- c(rep(0, s), base.X[i,])
      f.con2 <- cbind(base.Y, -base.X)
      f.con <- rbind(f.con1, f.con2, ar.l)
      f.dir <- c("==", rep("<=", nf), ar.dir)
      f.rhs <- c(1, rep(0, nf), ar.r)
      if(rts == 2){
        f.obj <- c(f.obj, -1, 1)
        z.tmp <- c(0, rep(1, nf), rep(0, nrow(ar.l)))
        z <- cbind(-z.tmp, z.tmp)
        f.con <- cbind(f.con, z)
      }
      if(dual == FALSE){
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
      else{
        re.tmp <- lp("max", f.obj, f.con, f.dir, f.rhs, compute.sens =
                     TRUE)$duals[1:(1+nf+m+s)]
        ## 여기에 q1, q1를 더하고 빼야 하는 작업을 해야 함.
        
      }
    }

    if(orientation == 2){
      f.obj <- c(rep(0, s), rep(base.X[i,], m))
      f.con1 <- c(base.Y[i,], rep(0, m))
      f.con2 <- cbind(base.Y, -base.X)
      f.con <- rbind(f.con1, f.con2, ar.l)
      f.dir <- c("==", rep("<=", nf), ar.dir)
      f.rhs <- c(1, rep(0, nf), ar.r)
      if(rts == 2){
        f.obj <- c(f.obj, -1, 1)
        z.tmp <- c(0, rep(1, nf), rep(0, nrow(ar.l)))
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


