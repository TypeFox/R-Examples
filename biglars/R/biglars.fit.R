biglars.fit <- 
function(x, y, type = "lasso", removeColumns = TRUE, eps = sqrt(.Machine$
  double.eps), blockSize = NULL, maxStages = NULL)
{
  #
  # This function calls a fitting routine, for "lar", "lasso" or
  # "stepwise" regression.
  #
  # x is a numeric matrix which does not include an explicit intercept
  # y is a numeric response
  #
  itype <- charmatch(type, c("lasso", "lar", "stepwise"))
  vecnorm <- function(x)
    sqrt(sum(x * x))

  qrFac <- qrBlockApply(x, y, blockSize)

  nrowx <- nrow(x)
  ncolx <- ncols <- ncol(x)
  colNames <- colnames(x)
  if(any(duplicated(colNames)))
    stop("column names must be unique")
  if(!length(colNames)) {
    colNames <- as.character(1:ncolx)
  }
  FUN <- switch(type,
    lasso = "biglars.fit.lasso",
    lar = "biglars.fit.lar",
    stepwise = "biglars.fit.stepwise",
    stop("unrecognized type"))
  out <- do.call(FUN, c(list(R = qrFac$R[-1, -1], Qty = qrFac$Qty[-1], 
    removeColumns = removeColumns, eps = eps, maxStages = maxStages)))
  coef <- rbind(0, out$coef)

  # INTERCEPT <- apply(coef, 1, function(b, x, y)
  # sum(y - x %*% b), x = x, y = y)/nrowx
  INTERCEPT <- apply(coef, 1, function(u, alpha, v)
  alpha - sum(u * v), alpha = qrFac$Qty[1], v = qrFac$R[1, -1]) * (qrFac$R[1, 1
    ]/nrowx)
  coef <- cbind(INTERCEPT, coef)
  colnames(coef) <- c("(Intercept)", colNames)
  RSS <- apply(coef, 1, function(b, Qty, R)
  vecnorm(Qty - R %*% b)^2, Qty = qrFac$Qty, R = qrFac$R)
  rownames(out$moves) <- colnames(x)[abs(out$moves[, "Var"])]
  structure(list(coefficients = coef, moves = out$moves, RSS = RSS), class = 
    "biglars")
}

