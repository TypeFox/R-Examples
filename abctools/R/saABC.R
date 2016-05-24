saABC <- function(theta, X, plot=TRUE) {
  if(nrow(theta)!=nrow(X)) stop("Number of rows in theta and X should be equal")
  X <- as.matrix(X)
  d <- ncol(theta)
  p <- ncol(X)
  n <- nrow(X)
  if(n<=p) stop("Too few simulated datasets to perform linear regression")
  if (is.null(colnames(theta))) {
    parnames <- paste("Parameter", 1:d)
  } else {
    parnames <- colnames(theta)
  }
  if (is.null(colnames(X))) {
    varnames <- paste("X", 1:p, sep="")
  } else {
    varnames <- colnames(X)
  }
  X <- cbind(1,X) ##Adds column for intercept

  if (plot) par(mfrow=c(1,d))
  B0 <- c()
  B <- matrix(nrow=d, ncol=p)
  BICs <- c()
  for (i in 1:d) { 
    reg <- lm.fit(X,theta[,i])
    class(reg) <- "lm" ##Needed to allow BIC calculation below
    if (plot) plot(theta[,i], reg$fitted.values,
                   xlab="True value", ylab="Fitted value",
                   main=parnames[i], pch=16) 
    B0[i] <- reg$coefficients[1]
    B[i,] <- reg$coefficients[-1]
    BICs[i] <- AIC(reg, k=log(n))
  }
  if (any(is.na(B))) warning("Linear regression problems: ill-conditioning?")
  names(B0) <- parnames
  names(BICs) <- parnames
  rownames(B) <- parnames
  colnames(B) <- varnames
  return(list(B0=B0,B=B,BICs=BICs))
}
