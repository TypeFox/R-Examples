`gini` <- function (var)
 {
  n <- length(var)
  vsort <- sort(var)

  Xk <- vsort

  f <- rep(1,length(Xk))
  f <- f/n

  F <- rep(0,length(Xk),1)

  for (i in 1:length(Xk))
   {
    F[i] <- length(which(var <= Xk[i]))
   }

  F <- F/n

  g <- (Xk*f)/(sum(var)/n)
  G <- cumsum(g)

  F1 <- matrix(rep(t(f),length(f)),ncol=dim(t(f))[2],byrow=FALSE)
  F2 <- matrix(rep(t(f),length(f)),ncol=dim(t(f))[2],byrow=TRUE)
  X1 <- matrix(rep(t(Xk),length(Xk)),ncol=dim(t(Xk))[2],byrow=FALSE)
  X2 <- matrix(rep(t(Xk),length(Xk)),ncol=dim(t(Xk))[2],byrow=TRUE)

  T <- (F1*F2)*abs(X1-X2)
  gini <- sum(T)/(2*sum(var)/n)

  F <- c(0,F)
  G <- c(0,G)

return(list(f=f, F=F, g=g, G=G, gini=gini))

}

