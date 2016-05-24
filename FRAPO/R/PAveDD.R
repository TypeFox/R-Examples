##
## Function for optimising a constrained average draw down portfolio
##
PAveDD <- function(PriceData, AveDD = 0.1, softBudget = FALSE, ...){
  if(is.null(dim(PriceData))){
    stop("Argument for 'PriceData' must be rectangular.\n")
  }
  if(any(is.na(PriceData))){
    stop("NA-values contained in object for 'PriceData'.\n")
  }
  if(AveDD <= 0 || AveDD >= 1){
    stop("Argument for 'AveDD' must be in the interval (0, 1).\n")
  }
  call <- match.call()
  RC <- as.matrix(returnseries(PriceData, method = "discrete", percentage = FALSE, compound = TRUE))
  rownames(RC) <- NULL
  N <- ncol(RC)
  J <- nrow(RC)
  w <- rep(0, N) ## weights
  u <- rep(0, J) ## high-watermark
  v <- rep(0, J) ## draw downs
  x <- c(w, u, v)
  ## Defining objective (end-wealth)
  obj <- c(as.numeric(RC[J, ]), rep(0, J), rep(0, J))
  ## a1: constraint that weights are positive
  a1 <- cbind(diag(N), matrix(0, nrow = N, ncol = 2 * J))
  d1 <- rep(">=", N)
  b1 <- rep(0, N)
  ## a2: budget constraint
  a2 <- c(rep(1, N), rep(0, 2 * J))
  ifelse(softBudget, d2 <- "<=", d2 <- "==")
  b2 <- 1
  ## a3: draw-down constraint (1) assigning summands to v
  a3 <- cbind(-1 * RC, diag(J), -1 * diag(J))
  d3 <- rep("==", J)
  b3 <- rep(0, J)
  ## a4: defining average constraint
  a4 <- c(rep(0, N), rep(0, J), rep(1 / J, J))
  d4 <- "<="
  b4 <- AveDD
  ## a5: draw-down constraint (2)
  a5 <- cbind(-1 * RC, diag(J), matrix(0, nrow = J, ncol = J))
  d5 <- rep(">=", J)
  b5 <- rep(0, J)
  ## a6: draw-down constraint (3)
  D1 <- -1.0 * diag(J)
  udiag <- embed(1:J, 2)[, c(2, 1)] 
  D1[udiag] <- 1
  a6 <- cbind(matrix(0, ncol = N, nrow = J), D1, matrix(0, ncol = J, nrow = J))
  a6 <- a6[-J, ]
  d6 <- rep(">=", J-1)
  b6 <- rep(0, J-1)
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
  ## Creating object PortMdd, inherits from PortSol
  weights <- opt$solution[1:N]
  names(weights) <- colnames(PriceData)
  dd <- timeSeries(opt$solution[(N + J + 1):(N + J + J)], charvec = rownames(PriceData))
  obj <- new("PortAdd", weights = weights, opt = opt, type = "average draw-down", call = call, AveDD = mean(dd), DrawDown = dd)
  return(obj)
}
