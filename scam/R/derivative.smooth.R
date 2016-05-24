
## function to get the derivative of the univariate constrained smooth term...


derivative.smooth <- function(object,smooth.number,x){
## object - a fitted scam object
## smooth.number - number of smooth term which derivative has to be found
## x - covariate values

  sn <- smooth.number
  m <- object$smooth[[sn]]$m
  q <- object$smooth[[sn]]$bs.dim
  n <- length(x)
     xk<-rep(0,q+m+2)
     xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
     for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
     for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  h <- (max(x)-min(x))/(q-m-1)
  Sig <- object$smooth[[sn]]$Sigma 

  first <- object$smooth[[sn]]$first.para
  last <- object$smooth[[sn]]$last.para
  beta.t <- object$coefficients.t[first:last]
  gamma <- Sig%*%beta.t
  delta.gamma <- diff(gamma, differences=1)
#  get model matrix for derivative function, of one order less -------------

  X1 <- splineDesign(xk,x,ord=m+1) # ord is by one less for the derivative
  deriv <- (X1[,2:(q-1)]%*%delta.gamma)/h ## derivative function for the monotone smooth
  
  ## calculating standard erors for getting CI...

  Sig1 <- matrix(0,q,q) 
  Sig1[,1] <- rep(1,q)
  Sig1[2:q,2:q] <- Sig 
  P <- diff(diag(q),difference=1)
  X <- (X1[,1:(q-1)]%*%P%*%Sig1)/h
  X <- sweep(X,2,colMeans(X))
  Vp <- object$Vp.t[c(1,first:last),c(1,first:last)] 
  Vp.c <- Vp
  Vp.c[,1] <- rep(0,nrow(Vp))
  Vp.c[1,] <- rep(0,ncol(Vp))
  se.fit <- sqrt(rowSums((X%*%Vp.c)*X))  

### trying another variant...

  P <- diff(diag(q-1),difference=1)
  X <- (X1[,1:(q-2)]%*%P%*%Sig)/h
  X <- sweep(X,2,colMeans(X))
  Vp <- object$Vp.t[first:last,first:last] 
  se.fit1 <- sqrt(rowSums((X%*%Vp)*X))  

list(deriv=deriv,se.deriv=se.fit,se.deriv1=se.fit1)
}



