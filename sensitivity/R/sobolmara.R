# Sobol' indices estimation by matrix random permutation (Mara 2008)
# Bootstrap is actually not possible
#
# Bertrand Iooss 2013


sobolmara <- function(model = NULL, X1, ...) {
  p <- ncol(X1)
  n <- nrow(X1)

  XX <- matrix(1:n,ncol=p,nrow=n)
  RP <- apply(XX,2,sample)
  X2 <- X1
  for (j in 1:p) X2[,j] <- X1[RP[,j],j]

  X <- rbind(X1, X2) 
  
  x <- list(model = model, X1 = X1, RP = RP, X = X, call = match.call())
  class(x) <- "sobolmara"
  
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }
  
  return(x)
}

estim.sobolmara <- function(data, i = 1 : nrow(data), RP) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  p <- ncol(RP)
  V <- var(d[, 1])
  m2 <- mean(d[,1])^2
  VCE <- NULL
  
  for (j in 1:p){
    hoy <- 0
    hoy <- sum(d[RP[,j],1]*d[,2])
    VCE <- cbind(VCE, hoy / (n - 1) - m2)
  }
  c(V, VCE)
}


tell.sobolmara <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  p <- ncol(x$X1)
  n <- nrow(x$X1)
  
  data <- matrix(x$y, nrow = n)
  
  # estimation of the variances of the conditional expectations (V)
  
  V <- data.frame(original = estim.sobolmara(data,RP=x$RP))
  rownames(V) <- c("global", colnames(x$X1))
  
  # estimation of the Sobol' indices
  
  S <- V[2:(p + 1), 1, drop = FALSE] / V[1,1] 
  rownames(S) <- colnames(x$X1)
  
  # return
  x$V <- V
  x$S <- S
  
  for (i in return.var) {
    x[[i]] <- get(i)
  }
  
  assign(id, x, parent.frame())
}


print.sobolmara <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      cat("\nSobol indices\n")
      print(x$S)
    }
  } else {
    cat("(empty)\n")
  }
}

plot.sobolmara <- function(x, ylim = c(0, 1), ...) {
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
  }
}
