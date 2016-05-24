PMinCDaR <- function(PriceData, alpha = 0.95, softBudget = FALSE, ...){
  if(is.null(dim(PriceData))){
    stop("Argument for 'PriceData' must be rectangular.\n")
  }
  if(any(is.na(PriceData))){
    stop("NA-values contained in object for 'PriceData'.\n")
  }
  if(alpha <= 0 || alpha >= 1){
    stop("Argument for 'alpha' must be in the interval (0, 1).\n")
  }
  call <- match.call()
  RC <- as.matrix(returnseries(PriceData, method = "discrete", percentage = FALSE, compound = TRUE))
  rownames(RC) <- NULL
  N <- ncol(RC)
  J <- nrow(RC)
  w <- rep(0, N) ## weights
  u <- rep(0, J) ## high-watermark
  v <- rep(0, J) ## draw downs (above threshold)
  z <- 0         ## thresh
  x <- c(w, u, v, z)
  ## Defining objective (end-wealth)
  obj <- c(rep(0, N), rep(0, J), (1/J)*(1/(1-alpha)) * rep(1, J), 1)
  ## a1: constraint that weights are positive
  a1 <- cbind(diag(N), matrix(0, nrow = N, ncol = 2 * J + 1))
  d1 <- rep(">=", N)
  b1 <- rep(0, N)
  ## a2: budget constraint
  a2 <- c(rep(1, N), rep(0, 2 * J + 1))
  ifelse(softBudget, d2 <- "<=", d2 <- "==")
  b2 <- 1
  ## a3: draw-down constraint (1) assigning summands to v
  a3 <- cbind(RC, -1 * diag(J), diag(J), matrix(1, ncol = 1, nrow = J))
  d3 <- rep(">=", J)
  b3 <- rep(0, J)
  ## a4: Unequality such draw downs above thresh are positive
  a4 <- cbind(matrix(0, ncol = N, nrow = J),
              matrix(0, ncol = J, nrow = J),
              diag(J),
              matrix(0, ncol = 1, nrow = J))
  d4 <- rep(">=", J)
  b4 <- rep(0, J)
  ## a5: draw-down constraint (2)
  a5 <- cbind(-1 * RC, diag(J), matrix(0, nrow = J, ncol = J + 1))
  d5 <- rep(">=", J)
  b5 <- rep(0, J)
  ## a6: draw-down constraint (3)
  D1 <- -1.0 * diag(J)
  udiag <- embed(1:J, 2)[, c(2, 1)] 
  D1[udiag] <- 1
  a6 <- cbind(matrix(0, ncol = N, nrow = J),
              D1,
              matrix(0, nrow = J, ncol = J + 1))
  a6 <- a6[-J, ]
  d6 <- rep(">=", J - 1)
  b6 <- rep(0, J - 1)
  ## a7: draw-down constraint (4)
  a7 <- c(rep(0, N), 1, rep(0, J-1), rep(0, J), 0)
  d7 <- "=="
  b7 <- 0
  ## Combining restrictions
  Amat <- rbind(a1, a2, a3, a4, a5, a6, a7)
  Dvec <- c(d1, d2, d3, d4, d5, d6, d7)
  Bvec <- c(b1, b2, b3, b4, b5, b6, b7)
  ## Solving LP
  opt <- Rglpk_solve_LP(obj = obj, mat = Amat, dir = Dvec, rhs = Bvec, ...)
  if(opt$status != 0){
    warning(paste("GLPK had exit status:", opt$status))
  }
  ## Creating object PortMdd, inherits from PortSol
  weights <- opt$solution[1:N]
  names(weights) <- colnames(PriceData)
  equity <- matrix(apply(RC, 1, function(x) sum(x * weights)), ncol = 1)
  rownames(equity) <- rownames(PriceData)
  uvals <- opt$solution[(N + 1):(N + J)]
  dd <- as.timeSeries(uvals - equity)
  colnames(dd) <- "DrawDowns"
  z <- opt$solution[N + J + J + 1]
  CDaR <- mean(dd[dd >= z, 1])
  obj <- new("PortCdd", weights = weights, opt = opt, type = "Conditional Draw Down at Risk", call = call, CDaR = CDaR, thresh = z, DrawDown = dd)
  return(obj)
}
