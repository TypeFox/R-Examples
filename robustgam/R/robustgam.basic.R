# smooth.construct.tr.spec is copied from mgcv package
smooth.construct.tr.smooth.spec<-function(object,data,knots){
  # from mgcv package
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  nk<-object$bs.dim-m-1 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  x.shift <- mean(x) # shift used to enhance stability
  k <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(k)) {
    n<-length(x)
    k<-quantile(x[2:(n-1)],seq(0,1,length=nk+2))[2:(nk+1)]
  }
  if (length(k)!=nk) # right number of knots?
    stop(paste("there should be ",nk," supplied knots"))
  x <- x - x.shift # basis stabilizing shift
  k <- k - x.shift # knots treated the same!
  X<-matrix(0,length(x),object$bs.dim)
  for (i in 1:(m+1)) X[,i] <- x^(i-1)
  for (i in 1:nk) X[,i+m+1]<-(x-k[i])^m*as.numeric(x>k[i])
  object$X<-X # the finished model matrix
  if (!object$fixed) {
    object$S[[1]]<-diag(c(rep(0,m+1),rep(1,nk)))
  }
  object$rank<-nk  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  ## store "tr" specific stuff ...
  object$knots<-k;object$m<-m;object$x.shift <- x.shift
  ## get the centering constraint ...
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tr.smooth"  # Give object a class
  object
}

# Predict.matrix.tr.smooth is copied from mgcv package
Predict.matrix.tr.smooth<-function(object,data){
  # from mgcv package
  x <- data[[object$term]]
  x <- x - object$x.shift # stabilizing shift
  m <- object$m;     # spline order (3=cubic)
  k<-object$knots    # knot locations
  nk<-length(k)      # number of knots
  X<-matrix(0,length(x),object$bs.dim)
  for (i in 1:(m+1)) X[,i] <- x^(i-1)
  for (i in 1:nk) X[,i+m+1] <- (x-k[i])^m*as.numeric(x>k[i])
  X # return the prediction matrix
}

# Huber loss function
Huber <- function(x,c){
  return( (abs(x)<=c)*x^2/2 + (abs(x)>c)*(c*abs(x)-c^2) )
}

Huber.deriv <- function(x,c){
  return( (abs(x)<=c)*x + (abs(x)>c)*(c*sign(x)) )
}

# robust estimation of sigma
robust.sigma <- function(res){
  #using the median absolute deviation of these residuals
  return(median(abs(res-median(res))))
}

# for gam fitting
mat.sqrt<- function(S){
  d<-eigen(S,symmetric=TRUE)
  if (sum(d$values<0)) {
  	index <- d$values<0
  	#cat("mat.sqrt: there are negative values:",d$values[index],"\n")
  	d$values[index]<- 0
  }
  rS<-d$vectors%*%diag(d$values^0.5)%*%t(d$vectors)
}

# Fisher consistency correction
## Poisson case
expect.poisson <- function(m, c, sqrtVar,...){
  k1 <- m-c*sqrtVar
  k2 <- m+c*sqrtVar
  return(sqrtVar*(dpois(ceiling(k1-1),m)-dpois(floor(k2),m))+c*(1-ppois(k2,m)-ppois(ceiling(k1-1),m)))
}

## Binomial (Bernoulli) case
expect.binomial <- function(m, c, sqrtVar,...){
  return( (Huber.deriv((1-m)/sqrtVar, c)*m+Huber.deriv(-m/sqrtVar,c)*(1-m)) )
}


expect.gaussian <- function(m, c, sqrtVar, ...){
  return(rep(0,length(m)))
}
