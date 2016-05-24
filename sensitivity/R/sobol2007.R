# Sobol' indices estimation (Sobol 2007 - Saltelli 2010)
#
# Author: Bertrand Iooss 2012


sobol2007 <- function(model = NULL, X1, X2, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)))
    stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1, X2)
  for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, nboot = nboot, conf = conf, X = X,
            call = match.call())
  class(x) <- "sobol2007"
  
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }
  
  return(x)
}


estim.sobol2007 <- function(data, i = 1 : nrow(data)) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  V <- var(d[, 1])
  VCE <- (colSums((d[, - c(1, 2)] - d[,1]) * d[, 2]) / (n - 1))
  VCE.compl <- (colSums((d[,1] - d[, - c(1, 2)]) * d[,1]) / (n - 1))
  c(V, VCE, VCE.compl)
}


tell.sobol2007 <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))

  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }

  p <- ncol(x$X1)
  n <- nrow(x$X1)

  data <- matrix(x$y, nrow = n)

  # estimation of the partial variances (V, D1 and Dt)
  
  if (x$nboot == 0){
    V <- data.frame(original = estim.sobol2007(data))
  }
  else{
    V.boot <- boot(data, estim.sobol2007, R = x$nboot)
    V <- bootstats(V.boot, x$conf, "basic")
  }
  rownames(V) <- c("global", colnames(x$X1), paste("-", colnames(x$X1), sep = ""))

  # estimation of the Sobol' indices (S1 and St)

  if (x$nboot == 0) {
    S <- V[2:(p + 1), 1, drop = FALSE] / V[1,1]
    T <- V[(p + 2):(2 * p + 1), 1, drop = FALSE] / V[1,1]
  } else {
    S.boot <- V.boot
    S.boot$t0 <- V.boot$t0[2:(p + 1)] / V.boot$t0[1]
    S.boot$t <- V.boot$t[,2:(p + 1)] / V.boot$t[,1]
    S <- bootstats(S.boot, x$conf, "basic")
    
    T.boot <- V.boot
    T.boot$t0 <- V.boot$t0[(p + 2):(2 * p + 1)] / V.boot$t0[1]
    T.boot$t <- V.boot$t[,(p + 2):(2 * p + 1)] / V.boot$t[,1]
    T <- bootstats(T.boot, x$conf, "basic")
  }
  rownames(S) <- colnames(x$X1)
  rownames(T) <- colnames(x$X1)

  # return
  x$V <- V
  x$S <- S
  x$T <- T

  for (i in return.var) {
    x[[i]] <- get(i)
  }

  assign(id, x, parent.frame())
}


print.sobol2007 <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nTotal indices:\n")
    print(x$T)
  }
}


plot.sobol2007 <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}
