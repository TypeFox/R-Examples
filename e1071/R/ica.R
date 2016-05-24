ica <- function(X, lrate, epochs=100, ncomp=dim(X)[2], 
                      fun="negative")
  {
    if (!is.matrix(X))
      {
        if (is.data.frame(X))
          X <- as.matrix(X)
        else
          stop("ica: X must be a matrix or a data frame")
      }
    if (!is.numeric(X))
      stop("ica: X contains non numeric elements")
            
    m <- dim(X)[1]
    n <- dim(X)[2]

    Winit <- matrix(rnorm(n*ncomp), ncomp, n)
    W <- Winit

    if (!is.function(fun))
      {
        funlist <- c("negative kurtosis", "positive kurtosis",
                     "4th moment")
        p <- pmatch(fun, funlist)
        if (is.na(p))
          stop("ica: invalid fun")
        funname <- funlist[p]
        if (p == 1) fun <- tanh
        else if (p == 2) fun <- function(x) {x - tanh(x)}
        else if (p == 3) fun <- function(x) {sign(x)*x^2}
      }
    else funname <- as.character(substitute(fun))
    
    for (i in 1:epochs)
      for (j in 1:m)
        {
          x <- X[j,, drop=FALSE]
          y <- W%*%t(x)
          gy <- fun(y)
          W <- W + lrate*gy%*%(x-t(gy)%*%W)
        }
    colnames(W) <- NULL
    pr <- X%*%t(W)
    retval <- list(weights = W, projection = pr, epochs = epochs,
                fun = funname, lrate = lrate, initweights = Winit)
    class(retval) <- "ica"
    return(retval)
  }


print.ica <- function(x, ...)
  {
    cat(x$epochs, "Trainingssteps with a learning rate of", x$lrate, "\n")
    cat("Function used:", x$fun,"\n\n")
    cat("Weightmatrix\n")
    print(x$weights, ...)
  }

plot.ica <- function(x, ...) pairs(x$pr, ...)

