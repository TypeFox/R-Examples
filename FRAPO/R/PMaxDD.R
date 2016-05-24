##
## Function for optimising a constrained maximum draw down portfolio
##
PMaxDD <- function(PriceData, MaxDD = 0.1, softBudget = FALSE, ...){
  if(is.null(dim(PriceData))){
    stop("Argument for 'PriceData' must be rectangular.\n")
  }
  if(any(is.na(PriceData))){
    stop("NA-values contained in object for 'PriceData'.\n")
  }
  if(MaxDD <= 0 || MaxDD >= 1){
    stop("Argument for 'MaxDD' must be in the interval (0, 1).\n")
  }
  call <- match.call()
  RC <- as.matrix(returnseries(PriceData, method = "discrete", percentage = FALSE, compound = TRUE))
  rownames(RC) <- NULL
  N <- ncol(RC)
  J <- nrow(RC)
  w <- rep(0, N) ## weights
  u <- rep(0, J) ## high-watermark
  x <- c(w, u)
  ## Defining objective (end-wealth)
  obj <- c(as.numeric(RC[J, ]), rep(0, J))
  ## a1: constraint that weights are positive
  a1 <- cbind(diag(N), matrix(0, nrow = N, ncol = J))
  d1 <- rep(">=", N)
  b1 <- rep(0, N)
  ## a2: budget constraint
  a2 <- c(rep(1, N), rep(0, J))
  ifelse(softBudget, d2 <- "<=", d2 <- "==")
  b2 <- 1
  ## a3: draw-down constraint (1)
  a3 <- cbind(-1 * RC, diag(J))
  d3 <- rep("<=", J)
  b3 <- rep(MaxDD, J)
  ## a4: draw-down constraint (2)
  a4 <- a3
  d4 <- rep(">=", J)
  b4 <- rep(0, J)
  ## a5: draw-down constraint (3)
  D1 <- -1.0 * diag(J)
  udiag <- embed(1:J, 2)[, c(2, 1)] 
  D1[udiag] <- 1
  a5 <- cbind(matrix(0, ncol = N, nrow = J), D1)
  a5 <- a5[-J, ]
  d5 <- rep(">=", J-1)
  b5 <- rep(0, J-1)
  ## a6: draw-down constraint (4)
  a6 <- c(rep(0, N), 1, rep(0, J - 1))
  d6 <- "=="
  b6 <- 0  
  ## Combining restrictions
  Amat <- rbind(a1, a2, a3, a4, a5, a6)
  Dvec <- c(d1, d2, d3, d4, d5, d6)
  Bvec <- c(b1, b2, b3, b4, b5, b6)
  ## Solving LP
  opt <- Rglpk_solve_LP(obj = obj, mat = Amat, dir = Dvec, rhs = Bvec,
                        max = TRUE, ...)
  if(opt$status != 0){
    warning(paste("GLPK had exit status:", opt$status))
  }
  ## Creating object PortMaxDD, inherits from PortSol
  weights <- opt$solution[1:N]
  names(weights) <- colnames(PriceData)
  equity <- matrix(apply(RC, 1, function(x) sum(x * weights)), ncol = 1)
  rownames(equity) <- rownames(PriceData)
  uvals <- opt$solution[(N + 1):(N + J)]
  dd <- as.timeSeries(uvals - equity)
  colnames(dd) <- "DrawDowns"
  obj <- new("PortMdd", weights = weights, opt = opt, type = "maximum draw-down", call = call, MaxDD = max(dd), DrawDown = dd)
  return(obj)
}
