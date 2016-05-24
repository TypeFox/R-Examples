revenue.dea <- function(base = NULL, frontier = NULL, 
                        noutput = 1, output.price = NULL) {

  ## output.price: c(p1', p2', ..., ps')
  
  require(lpSolve)

  if(is.null(frontier))
    frontier <- base

  if(!is.null(base) & !is.null(frontier)){
    base <- as.matrix(base)
    frontier <- as.matrix(frontier)
  }

  if(ncol(base) != ncol(frontier))
    stop("Number of columns in base matrix and frontier matrix should be the same!")

  if(!is.vector(output.price))
    stop("Fixed output price (vector) must be provided by user!!")
  
  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)
  nf <- nrow(frontier)

  front.Y <- t(frontier[, 1:s])
  front.X <- t(frontier[, (s+1):(s+m)])
  base.Y <- t(base[, 1:s])
  base.X <- t(base[, (s+1):(s+m)])

  ps <- output.price
  each.revenue <- base.Y * output.price
  total.revenue <- apply(each.revenue, 2, sum)

  re <- data.frame(matrix(0, nrow = n, ncol = s + 6))
  names(re) <- c(paste("y.star", 1:s, sep = ""),
                 "OE.revenue", "TE.revenue", "Revealed.revenue", "OE", 
                 "AE", "TE")

  thetas <- dea(base = base, frontier = frontier, noutput = noutput,
                orientation = 2, rts = 1)[, 1]

  for(i in 1:n){
    f.obj <- c(ps, rep(0, nf))
    f.dir <- c(rep(">=", s), rep("<=", m))
    f.rhs <- c(rep(0, s), base.X[,i])

    f.con1 <- cbind(-diag(s), front.Y)
    f.con2 <- cbind(matrix(0, m, s), front.X)
    f.con <- rbind(f.con1, f.con2)

    re.tmp <- lp("max", f.obj, f.con, f.dir, f.rhs)

    max.revenue<- re.tmp$objval
    y.star <- re.tmp$solution[1:s]

    re[i, 1:s] <- y.star
    re[i, s+1] <- max.revenue
    re[i, s+2] <- thetas[i] * sum(base.Y[,i]*ps) # technical efficiency
    re[i, s+3] <- sum(base.Y[,i]*ps)
    re[i, s+4] <- re[i, s+3]/re[i, s+1]
    re[i, s+5] <- re[i, s+2]/re[i, s+1]
    re[i, s+6] <- 1/thetas[i]
  }

  return(re)
}
