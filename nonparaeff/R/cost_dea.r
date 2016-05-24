cost.dea <- function(base = NULL, frontier = NULL,
                     noutput = 1, input.price = NULL) {

  ## input.price: c(w1', w2', ..., wm')
  
  require(lpSolve)

  if(is.null(frontier))
    frontier <- base

  if(!is.null(base) & !is.null(frontier)){
    base <- as.matrix(base)
    frontier <- as.matrix(frontier)
  }

  if(ncol(base) != ncol(frontier))
    stop("Number of columns in base matrix and frontier matrix should be the same!")

  if(!is.vector(input.price))
    stop("Fixed input price (vector) must be provided by user!!")
  
  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)
  nf <- nrow(frontier)

  front.Y <- t(frontier[, 1:s])
  front.X <- t(frontier[, (s+1):(s+m)])
  base.Y <- t(base[, 1:s])
  base.X <- t(base[, (s+1):(s+m)])

  wm <- input.price
  each.cost <- base.X * input.price
  total.cost <- apply(each.cost, 2, sum)

  re <- data.frame(matrix(0, nrow = n, ncol = m + 6))
  names(re) <- c(paste("x.star", 1:m, sep = ""),
                 "OE.cost",
                 "TE.cost", "Revealed.cost", "OE", "AE", "TE")

  thetas <- dea(base = base, frontier = frontier, noutput = noutput,
                orientation = 1, rts = 1)[, 1]

  for(i in 1:n){
    f.obj <- c(wm, rep(0, nf))
    f.dir <- c(rep("<=", m), rep(">=", s))
    f.rhs <- c(rep(0, m), base.Y[,i])

    f.con1 <- cbind(-diag(m), front.X)
    f.con2 <- cbind(matrix(0, s, m), front.Y)
    f.con <- rbind(f.con1, f.con2)

    re.tmp <- lp("min", f.obj, f.con, f.dir, f.rhs)

    min.cost <- re.tmp$objval
    x.star <- re.tmp$solution[1:m]

    re[i, 1:m] <- x.star
    re[i, m+1] <- min.cost # m + 1; OE cost
    re[i, m+2] <- thetas[i]*sum(base.X[,i]*wm) # m + 2; TE cost
    re[i, m+3] <- sum(base.X[,i]*wm) ## m+3; Revealed cost
    re[i, m+4] <- re[i, m+1]/re[i, m+3] ## m+4: OE
    re[i, m+5] <- re[i, m+4]/thetas[i]
    re[i, m+6] <- thetas[i] ## TE
  }

  return(re)
}
