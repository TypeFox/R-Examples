################################################################
##  Various dependence measures
################################################################

.covariance <- function(t,u){sum(t*u)}

.spearman <- function(x,y){sum(rank(x)*rank(y))}

.kendall <- function(t,u) {
  n <- length(t)
  sum(outer(t,t,FUN="-")*outer(u,u,FUN="-")>0)
}

corrob <- function(t,u) covRob(cbind(t,u), corr = TRUE, distance = FALSE, estim = "pairwiseGK")$cov[1,2]
# covRob() comes from package robust
covrob <- function(t, u) ((scaleTau2(t + u))^2 - (scaleTau2(t - u))^2)/4
# scaleTau2() comes from package robustbase

dcov <- function(x,y,Cpp=TRUE) {
# Distance covariance from Gabor J. Szekely et al., Annals of Stat, 2007, vol 35 (6), p.2769-2794
# Warning: Only valid to compute the distance covariance for two random variables X and Y
# This means that X and Y cannot be random Vectors.

  if (is.matrix(x)) if (ncol(x)>1) stop("Consider using the dcov() function in package energy.")
  if (is.matrix(y)) if (ncol(y)>1) stop("Consider using the dcov() function in package energy.")
  n <- length(x)
  if (length(y) != n) stop("x and y should have the same length")
  
  if (Cpp) {
    Vup <- .C("dcov",as.double(x),as.double(y),as.integer(n),Vup=as.double(0))$Vup
  } else {
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
    
    Vup <- sqrt(mean(A*B))
  }
  
  return(Vup)
  
}
