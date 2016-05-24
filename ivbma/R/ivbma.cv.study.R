predict.ivbma <- function(a,X,W)
  {
    n <- dim(X)[1]
    odens <- dim(a$rho)[1]
    V <- as.matrix(cbind(X,W))
    Y.hat <- matrix(0,100 * odens,1)
    for(i in 1:odens)
      {
        for(j in 1:n)
          {
            v <- a$Sigma[1,1,i]
            eps <- rnorm(100,0,sqrt(v))
            rho <- a$rho[i,]
            Y.mean <- V[j,] %*% rho
            l <- ( (i - 1) * 100 + 1):( (i - 1) * 100 + 100)
            Y.hat[l,j] <- Y.mean + eps
          }
      }
    return(Y.hat)
  }

ivbma.score <- function(b, y)
  {
    y.bar <- mean(b)
    SE <- (y - y.bar)^2
    AE <- abs(y - median(b))
    VAR <- var(b)
    CRPS <- NULL
    nn <- length(b)
    f <- mean( abs(y - b))
    f <- f - 1/2 * mean( abs(b[1:(nn - 1)] - b[2:nn]))
    CRPS <- f
    r <- c(SE,AE,VAR,CRPS)
    names(r) <- c("SE","AE","VAR","CRPS")
    return(r)
  }

ivbma.cv.study <- function(d,...)
  {
    n <- length(d$Y)
    R <- matrix(0,n,4)
    for(j in 1:n)
      {
        a <- ivbma(d$Y[-j],d$X[-j,],d$Z[-j,],d$W[-j,],...)
        b <- predict(a,d$X[j,],d$W[j,])
        R[j,] <- ivbma.score(b,d$Y[j])
      }
    return(R)
  }
