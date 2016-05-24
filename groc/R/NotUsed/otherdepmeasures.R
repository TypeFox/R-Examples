# Functions to compute Deheuvels measure of dependence
#Dn = function(s,t,n) { # Genest et. al page 276
#  (2*n+1)/(6*n)+s*(s-1)/(2*n*(n+1))+t*(t-1)/(2*n*(n+1))-max(s,t)/(n+1)
#}

# If we want to use outer in Bn() we must code Dn as below!!
Dn <- function(s,t,n) { # Genest et. al page 276
  res <- rep(NA,length(s))
  for (i in 1:length(s)) {
    res[i] <- (2*n+1)/(6*n)+s[i]*(s[i]-1)/(2*n*(n+1))+t[i]*(t[i]-1)/(2*n*(n+1))-max(s[i],t[i])/(n+1)
  }
  return(res)
}


deheuvels <- function(r,s,ranking=TRUE,Cpp=TRUE) { # Genest et. al page 276

  n <- length(r)
  if (length(s) != n) stop("r and s should have the same length")
  
  if (Cpp) {
    resBn <- 0
    out <- .C("deheuvels",as.double(r),as.double(s),as.integer(ranking),as.integer(n),resBn=as.double(resBn),PACKAGE="grolv")
    res <- out$resBn   
  } else {
    if (ranking) { # We should set ranking=FALSE only if r and s are already vectors of ranks
      r <- rank(r,ties.method ="first")
      s <- rank(s,ties.method ="first")
    }
    res <- sum(outer(r,r,Dn,n)*outer(s,s,Dn,n))/n
  }
  return(res)
}

correlation <- function(t,u){sum(t*u)/sqrt(sum(t^2)*sum(u^2))}

D.rho <- function(t,u) {
  n <- length(t)
  (12* spearman(t,u)/(n*(n^2-1))-3*(n+1)/(n-1))*mad(t)*mad(u)
}

D.kendall <- function(t,u) {n <- length(t);(2*kendall(t,u)/(n*(n-1))-1)*mad(t)*mad(u)}

dcor <- function(x,y) {
# Distance correlation from Gabor et al., Annals of Stat, 2007, vol 35 (6), p.2769-2794

  n <- length(x)
  if (length(y) != n) stop("x and y should have the same length")
  
  a <- outer(as.vector(x),as.vector(x),FUN=function(y,x) abs(y-x))
  mean.ak. <- rowMeans(a)
  mean.a.l <- colMeans(a)
  mean.a <- mean(a)
  A <- sweep(a,MARGIN=1,STATS=mean.ak.,FUN="-")
  A <- sweep(A,MARGIN=2,STATS=mean.a.l,FUN="-")
  A <- A + mean.a
  
  b <- outer(as.vector(y),as.vector(y),FUN=function(y,x) abs(y-x))
  mean.bk. <- rowMeans(b)
  mean.b.l <- colMeans(b)
  mean.b <- mean(b)
  B <- sweep(b,MARGIN=1,STATS=mean.bk.,FUN="-")
  B <- sweep(B,MARGIN=2,STATS=mean.b.l,FUN="-")
  B <- B + mean.b

  Vup <- sum(A*B)
  Vx <- sum(A^2)
  Vy <- sum(B^2)

  R2 <- Vup/sqrt(Vx*Vy)

  return(R2)
  
}

dcov.norm <- function(x,y) {
# Distance correlation from Gabor et al., Annals of Stat, 2007, vol 35 (6), p.2769-2794

  n <- length(x)
  if (length(y) != n) stop("x and y should have the same length")
  
  a <- outer(as.vector(x),as.vector(x),FUN=function(y,x) abs(y-x))
  mean.ak. <- rowMeans(a)
  mean.a.l <- colMeans(a)
  mean.a <- mean(a)
  A <- sweep(a,MARGIN=1,STATS=mean.ak.,FUN="-")
  A <- sweep(A,MARGIN=2,STATS=mean.a.l,FUN="-")
  A <- A + mean.a
  
  b <- outer(as.vector(y),as.vector(y),FUN=function(y,x) abs(y-x))
  mean.bk. <- rowMeans(b)
  mean.b.l <- colMeans(b)
  mean.b <- mean(b)
  B <- sweep(b,MARGIN=1,STATS=mean.bk.,FUN="-")
  B <- sweep(B,MARGIN=2,STATS=mean.b.l,FUN="-")
  B <- B + mean.b

  nS2 <- (sum(a)*sum(b))/(n^3)

  Vup <- sum(A*B)

  res <- Vup/nS2
  
  return(res)
  
}

D.Qn <- function(t,u) ((Qn(t+u))^2-(Qn(t-u))^2)/4
