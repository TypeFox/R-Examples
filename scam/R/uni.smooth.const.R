

#####################################################
### Adding Monotone increasing SCOP-spline construction 
######################################################


smooth.construct.mpi.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing smooth
{ 
   #require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ", nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   
     # elements of matrix Sigma for increasing smooth
  for (i in 1:(q-1))  Sig[i,1:i]<-1
  X <- X1[,2:q]%*%Sig # model submatrix for the monotone term
  object$X <- X # the finished model matrix
  object$P <- list()
  object$S <- list()
  object$Sigma <- Sig

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  b <- rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)

  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  class(object)<-"mpi.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpi` smooth class... 


Predict.matrix.mpi.smooth<-function(object,data)
## prediction method function for the `mpi' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   # Define Matrix Sigma
  # elements of matrix Sigma for increasing smooth
  for (i in 1:q)  Sig[i,1:i] <- 1
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


########################################################
### Adding Monotone decreasing SCOP-spline construction 
########################################################


smooth.construct.mpd.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubic splines
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n <- length(x)
       xk <- rep(0,q+m+2)
       xk[(m+2):(q+1)] <- seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i] <- xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i] <- xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
      # elements of matrix Sigma for decreasing smooth
  for (i in 1:(q-1))  Sig[i,1:i]<- -1
  X <- X1[,2:q]%*%Sig # model submatrix for the monotone term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-1),difference=1)
    object$P[[1]] <- P
    object$S[[1]] <- crossprod(P)
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mpd.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mpd` smooth class 

Predict.matrix.mpd.smooth<-function(object,data)
## prediction method function for the `mpd' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  # elements of matrix Sigma for decreasing smooth...
  Sig[,1] <- 1
  for (i in 2:q)  Sig[i,2:i]<- -1 
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}




##############################################################
### Smooth constructor for the mixed constrainted smooths ......
##############################################################

###############################################################
### Adding Monotone decreasing & concave P-spline construction ###############################################################


smooth.construct.mdcv.smooth.spec<- function(object, data, knots)
## construction of the monotone decreasing and concave smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default for cubis spline
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
     { n<-length(x)
       xk<-rep(0,q+m+2)
       xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
       for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
       for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
     }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & concave smooth
  for (i in 1:(q-1)) Sig[i:(q-1),i]<--c(1:(q-i))
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcv` smooth class.... 

Predict.matrix.mdcv.smooth<-function(object,data)
## prediction method function for the `mdcv' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q)  Sig[i:q,i] <- -c(1:(q-i+1))
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}



###############################################################
### Adding decreasing & convex SCOP-spline construction ################################################################


smooth.construct.mdcx.smooth.spec<- function(object, data, knots)
##  the constructor for the monotone decreasing and convex smooth
{ 
  #require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone decreasing & convex smooth
  Sig[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig[i,1:(q-i)] <- -i;
         Sig[i,(q-i+1):(q-1)] <- -c((i-1):1)
  }
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"mdcx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `mdcx` smooth class.... 


Predict.matrix.mdcx.smooth<-function(object,data)
## the prediction method for the `mdcx' smooth class
{ x <- data[[object$term]]
  m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- -rep(1,q-1)
  for (i in 2:(q-1)) {
         Sig1[i,1:(q-i)]<--i;
         Sig1[i,(q-i+1):(q-1)]<--c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


################################################################
### Adding monotone increasing & concave SCOP-spline construction ################################################################


smooth.construct.micv.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and concave smooth
{ 
 # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & concave smooth
  Sig[1,]<-rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig[i,1:(q-i)]<-i;
       Sig[i,(q-i+1):(q-1)]<-c((i-1):1)
  }
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micv` smooth class....
 

Predict.matrix.micv.smooth<-function(object,data)
## prediction method function for the `micv' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig1 <- matrix(0,(q-1),(q-1))   
  Sig1[1,] <- rep(1,q-1)
  for (i in 2:(q-1)) {
       Sig1[i,1:(q-i)] <- i;
       Sig1[i,(q-i+1):(q-1)] <- c((i-1):1)
  }
  Sig [2:q,2:q] <- Sig1

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}



###############################################################
### Adding monotone increasing & convex SCOP-spline construction ################################################################


smooth.construct.micx.smooth.spec<- function(object, data, knots)
## construction of the monotone increasing and convex smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m < 1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Matrix Sigma
     # for monotone increasing & convex smooth
  for (i in 1:(q-1)) Sig[i:(q-1),i]<-c(1:(q-i))
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
 
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"micx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `micx` smooth class...
 

Predict.matrix.micx.smooth<-function(object,data)
## prediction method function for the `micx` smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  for (i in 2:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


############################################################
### Smooth construct for the convex/concave smooths ......
###########################################################

###########################################################
### Adding concave SCOP-spline construction *************
###########################################################

smooth.construct.cv.smooth.spec<- function(object, data, knots)
## construction of the concave smooth
{ 
 # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))  # Define Sigma for concave smooth
  Sig[1:(q-1),1]<- c(1:(q-1))
  for (i in 2:(q-1)) Sig[i:(q-1),i]<--c(1:(q-i))
  
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"cv.smooth"  # Give object a class
  object
}


## Prediction matrix for the `cv` smooth class******************

Predict.matrix.cv.smooth<-function(object,data)
## prediction method function for the `cv' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- -c(1:(q-i+1))

  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}


###########################################################
### Adding convex SCOP-spline construction *************
##########################################################

smooth.construct.cx.smooth.spec<- function(object, data, knots)
## construction of the convex smooth
{ 
  # require(splines)
  m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default 
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## default
  q <- object$bs.dim 
  nk <- q+m+2 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]]  ## the data
  xk <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through data
  { n<-length(x)
    xk<-rep(0,q+m+2)
    xk[(m+2):(q+1)]<-seq(min(x),max(x),length=q-m)
    for (i in 1:(m+1)) {xk[i]<-xk[m+2]-(m+2-i)*(xk[m+3]-xk[m+2])}
    for (i in (q+2):(q+m+2)) {xk[i]<-xk[q+1]+(i-q-1)*(xk[m+3]-xk[m+2])}
  }
  if (length(xk)!=nk) # right number of knots?
  stop(paste("there should be ",nk," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m+2)
  # use matrix Sigma and remove the first column for the monotone smooth
  Sig <- matrix(0,(q-1),(q-1))   # Define Sigma for convex smooth
  Sig[1:(q-1),1]<- -c(1:(q-1))
  for (i in 2:(q-1)) Sig[i:(q-1),i]<- c(1:(q-i))
  
  X <- X1[,2:q]%*%Sig # model submatrix for the constrained term
  object$X<-X # the finished model matrix
  object$Sigma <- Sig
  object$P <- list()
  object$S <- list()

  if (!object$fixed) # create the penalty matrix
  { P <- diff(diag(q-2),difference=1)
    object$P[[1]] <- P
    S <- matrix(0,q-1,q-1)
    S[2:(q-1),2:(q-1)] <- crossprod(P)
    object$S[[1]] <- S
  }
  b<-rep(1,q-1) # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  
  ## get model matrix for 1st and 2nd derivatives of the smooth...
  h <- (max(x)-min(x))/(q-m-1) ## distance between two adjacent knots
  object$Xdf1 <- splineDesign(xk,x,ord=m+1)[,2:(q-1)]/h ## ord is by one less for the 1st derivative
  object$Xdf2 <- splineDesign(xk,x,ord=m)[,2:(q-2)]/h^2 ## ord is by two less for the 2nd derivative

  object$knots <- xk;
  object$m <- m;
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"cx.smooth"  # Give object a class
  object
}


## Prediction matrix for the `cx` smooth class *************************

Predict.matrix.cx.smooth<-function(object,data)
## prediction method function for the `cx' smooth class
{ m <- object$m+1; # spline order, m+1=3 default for cubic spline
  q <- object$df +1
  Sig <- matrix(0,q,q)   
  Sig[,1] <- rep(1,q)
  Sig[2:q,2]<- -c(1:(q-1))
  for (i in 3:q) Sig[i:q,i] <- c(1:(q-i+1))
  
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X <- spline.des(object$knots,x,m)$design
     X <- X%*%Sig 
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1))$design
     X <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X[ind,] <- spline.des(object$knots,x[ind],m)$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
     X <- X%*%Sig 
  }
  X
}





