#########################################################################
### Shape constrained smooth construct for bivariate terms......       ##
#########################################################################


####################################################################################### 
### Tensor product P-spline construction with double monotone decreasing constraint  ##
#######################################################################################

smooth.construct.tedmd.smooth.spec<- function(object, data, knots)
{ ## construction of the double monotone decreasing smooth surface
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  if (object$bs.dim[1]==-1)  # set the default values for q1 and q2
     {  q1 <- object$bs.dim[1] <- 7
        q2 <- object$bs.dim[2] <- 7
     }
  else if (length(object$bs.dim)==1)
       q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1))  q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2))  q2 <- object$bs.dim[2] <- 7
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
     {  n<-length(x)
        xk<-rep(0,q1+m[1]+2)
        xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
        for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
        for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
     }
  if (n != length(z))
       stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
     {  zk<-rep(0,q2+m[2]+2)
        zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
        for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
        for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
     }
  if (length(xk)!=nk1 || length(zk)!=nk2) # right number of knots?
       stop(paste("there should be ",nk1, " and ", nk2," supplied knots"))
  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
  for (j in 1:q2)  IS[j,1:j] <- -1
  IS1 <- matrix(0,q1,q1)   # Define submatrix of Sigma
  for (j in 1:q1)  IS1[j,1:j] <- 1
  Sig <- IS1%x%IS # Knonecker product to get Sigma
  Sig[,1] <- rep(1,ncol(Sig))
 
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-1),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2

  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]]) ## t(object$P[[2]])%*%object$P[[2]]

  b <- rep(1,q1*q2-1)  # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)

  ## store "tedmd" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tedmd.smooth"  # Give object a class
  object
}

####################################################################


Predict.matrix.tedmd.smooth <- function(object, data)
{ ## prediction method function for the `tedmd' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
  for (j in 1:q2)  IS[j,1:j] <- -1
  IS1 <- matrix(0,q1,q1)   # Define submatrix of Sigma
  for (j in 1:q1)  IS1[j,1:j] <- 1
  Sig <- IS1%x%IS # Knonecker product to get Sigma
  Sig[,1] <- rep(1,ncol(Sig))
  X <- X%*%Sig
  X # return the prediction matrix
}


########################################################################
## function used for predict method to get marginal model submatrices ## 
## with linear extrapolation if needed                                ##
########################################################################                                

marginal.linear.extrapolation <- function(object, data)
{ ## function to get marginal matrices used in predict method on bivariate SCOP-splines
  x <- data[[object$term[1]]]  
  z <- data[[object$term[2]]]  
  if (length(x) != length(z))
        stop ("arguments of smooth are not of the same dimension")
  m <- object$m + 1 ## vector of two components
  ## find spline basis inner knot range for 1st covariate, x...
  ll <- object$knots[[1]][m[1]+1];ul <- object$knots[[1]][length(object$knots[[1]])-m[1]]
  m[1] <- m[1] + 1
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X1 <- spline.des(object$knots[[1]],x,m[1])$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots[[1]],c(ll,ll,ul,ul),m[1],c(0,1,0,1))$design
     X1 <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X1[ind,] <- spline.des(object$knots[[1]],x[ind],m[1])$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- x < ll 
     if (sum(ind)>0) X1[ind,] <- cbind(1,x[ind]-ll)%*%D[1:2,]
     ind <- x > ul
     if (sum(ind)>0) X1[ind,] <- cbind(1,x[ind]-ul)%*%D[3:4,]
  }
  ## the same for 2nd maginal matrix...
  ## find spline basis inner knot range for 2nd covariate, z...
  ll <- object$knots[[2]][m[2]+1];ul <- object$knots[[2]][length(object$knots[[2]])-m[2]]
  m[2] <- m[2] + 1
  ind <- z<=ul & z>=ll ## data in range
  if (sum(ind)==n) { ## all in range
     X2 <- spline.des(object$knots[[2]],z,m[2])$design
  } else { ## some extrapolation needed 
     ## matrix mapping coefs to value and slope at end points...
     D <- spline.des(object$knots[[2]],c(ll,ll,ul,ul),m[2],c(0,1,0,1))$design
     X2 <- matrix(0,n,ncol(D)) ## full predict matrix
     if (sum(ind)> 0)  X2[ind,] <- spline.des(object$knots[[2]],z[ind],m[2])$design ## interior rows
     ## Now add rows for linear extrapolation...
     ind <- z < ll 
     if (sum(ind)>0) X2[ind,] <- cbind(1,z[ind]-ll)%*%D[1:2,]
     ind <- z > ul
     if (sum(ind)>0) X2[ind,] <- cbind(1,z[ind]-ul)%*%D[3:4,]
  }
  list(X1=X1, X2=X2)
}



####################################################################################### 
### Tensor product P-spline construction with double monotone increasing constraint  ##
#######################################################################################


smooth.construct.tedmi.smooth.spec <- function(object, data, knots)
{ ## construction of the double monotone increasing smooth surface
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  if (object$bs.dim[1]==-1)  # set the default values fro q1 and q2
     {  q1 <- object$bs.dim[1] <- 7
        q2 <- object$bs.dim[2] <- 7
     }
  else if (length(object$bs.dim)==1)
         q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
     {  n<-length(x)
        xk<-rep(0,q1+m[1]+2)
        xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
        for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
        for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
     }
  if (n != length(z))
       stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
      {  zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
      }
  if (length(xk)!=nk1 || length(zk)!=nk2) # right number of knots?
         stop(paste("there should be ",nk1, " and ", nk2," supplied knots"))

  #  get model matrix-------------
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  X2 <- splineDesign(zk,z,ord=m[2]+2)
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
      {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma
  IS2[1:q2,1] <- rep(1,q2)
    for (j in 2:q2)  IS2[j,2:j] <- 1
  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
  IS1[1:q1,1] <- rep(1,q1)
    for (j in 2:q1)  IS1[j,2:j] <- 1
  Sig <- IS1 %x% IS2 
  
  # apply identifiability constraint and get model matrix
  X <- X[,2:ncol(X)]%*%Sig[2:ncol(Sig),2:ncol(Sig)]  
  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  S <- list()
  I2<- diag(q2)
  I1 <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- I1
  S[[1]]<- Pm1%x%I2

  I2 <- diff(diag(q2-1),difference=1) 
  Pm2 <- matrix(0,q2-1,q2)
  Pm2[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%Pm2
  
  object$P <- list()
  object$P[[1]] <- S[[1]][2:nrow(S[[1]]),2:ncol(S[[1]])]
  object$P[[2]] <- S[[2]][2:nrow(S[[2]]),2:ncol(S[[2]])]
  object$S <- list()
  object$S[[1]] <- crossprod(object$P[[1]]) ## t(object$P[[1]])%*%object$P[[1]]
  object$S[[2]] <- crossprod(object$P[[2]])  ## t(object$P[[2]])%*%object$P[[2]]

  b <- rep(1,q1*q2-1)  # define vector of 0's & 1's for model parameters identification
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- diag(q1*q2-1) # identfiability constraint matrix
  object$Zc <- rbind(rep(0,ncol(object$Zc)),object$Zc)
  
  ## store "tedmi" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tedmi.smooth"  # Give object a class
  object
}


####################################################################


Predict.matrix.tedmi.smooth <- function(object, data)
{  ## prediction method function for the `tedmi' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for ( i in 1:n)
      {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
      }
  # get a matrix Sigma -----------------------
  IS2 <- matrix(0,q2,q2)   # Define marginal matrix of Sigma
  IS2[1:q2,1] <- rep(1,q2)
    for (j in 2:q2)  IS2[j,2:j] <- 1
  IS1 <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
  IS1[1:q1,1] <- rep(1,q1)
    for (j in 2:q1)  IS1[j,2:j] <- 1
  Sig <- IS1 %x% IS2 
  # get final model matrix
  X <- X %*%Sig
  X # return the prediction matrix
}



############################################################################ 
## Tensor product P-spline construction with single monotone decreasing   ##
## constraint wrt the first covariate                                     ##
############################################################################
                                                    
smooth.construct.tesmd1.smooth.spec<- function(object, data, knots)
{ ## construction of the double monotone increasing smooth surface
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      { m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
        object$p.order <- m
      }
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  object$p.order[is.na(object$p.order)] <- 2
  if (object$bs.dim[1]==-1)  # set the default values fro q1 and q2
     {  q1 <- object$bs.dim[1] <- 7
        q2 <- object$bs.dim[2] <- 7
     }
  else if (length(object$bs.dim)==1)
         {  q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
                   ## basis dimension is provided for both marginal smooths
            object$bs.dim <- rep(object$bs.dim, 2)
         }
  else {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
    
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")

  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
     {  n<-length(x)
        xk<-rep(0,q1+m[1]+2)
        xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
        for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
        for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
        knots[[object$term[1]]] <- xk
     }
  if (n != length(z))
       stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
     {   zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
     }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  #  get model matrix-------------
  # get marginal model matrices and penalties... 
  bm <- marginal.matrices.tesmi1.ps(x,z,xk,zk,m,q1,q2)
  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S      
  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
  IS[1:q1,1]<-1
  for (j in 2:q1)  IS[j,2:j] <- -1
  I <- diag(q2)
  Sig <- IS%x%I

  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D <- diag(q1*q2)
  D <- D[,-q2]
  D1 <- t(diff(diag(q2)))
  D[1:q2,1:(q2-1)] <- D1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D ##  t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D ## t(D)%*%S[[2]]%*%D

  b <- rep(1,q1*q2-1)  # define vector of 0's & 1's for model parameters identification
  b[1:(q2-1)] <- rep(0, q2-1) 
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmd1" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  if (is.null(zk))
     object$knots[[2]] <- rep(0,0,0)
  else object$knots[[2]] <- zk
  object$m <- m
      
  object$df <- ncol(object$X)     # maximum DoF (if unconstrained)
  class(object) <- "tesmd1.smooth"  # Give object a class
  object
}



###########################################################################

## Prediction matrix for the `tesmd1` smooth class *************************

Predict.matrix.tesmd1.smooth<-function(object,data)
{ ## prediction method function for the `tesmd1' smooth class
  if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
  IS[1:q1,1]<-1
  for (j in 2:q1)  IS[j,2:j] <- -1
  I <- diag(q2)
  Sig <- IS%x%I
  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}



############################################################################ 
## Tensor product P-spline construction with single monotone decreasing   ##
## constraint wrt the second covariate                                    ##
############################################################################


smooth.construct.tesmd2.smooth.spec<- function(object, data, knots)
## construction of the double monotone increasing smooth surface
{ 
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  
  if (length(object$p.order)==1)
      { m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
        object$p.order <- m
  }
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  object$p.order[is.na(object$p.order)] <- 2
  if (object$bs.dim[1]==-1) { # set the default values fro q1 and q2
      q1 <- object$bs.dim[1] <- 7
      q2 <- object$bs.dim[2] <- 7
  }
  else if (length(object$bs.dim)==1){
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
      object$bs.dim <- rep(object$bs.dim, 2)
  }
  else {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
  
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
         
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
       { zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
  }
  if (length(zk)!=nk2 ) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  if (length(xk)!=nk1) # right number of knots?
      stop(paste("there should be ",nk1," supplied knots for x"))
  
  #  get model matrix-------------
 # get marginal model matrices and penalties... 
  bm <- marginal.matrices.tesmi2.ps(x,z,xk,zk,m,q1,q2)
  
  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S  

  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
  
  IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
  IS[1:q2,1]<-1
  for (j in 2:q2)  IS[j,2:j] <- -1
  I <- diag(q1)
  Sig <- I%x%IS
    
  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D<- diag(q1*q2)
  D<-D[,-((q1-1)*q2+1)]
  ind <- rep(0,q1-1) # get index number for the cells to be changed to "-1"
  for (i in 1:(q1-1)){
        ind[i] <- (i-1)*q2+1
        D[ind[i],ind[i]] <- -1
  }
  for (i in 2:(q1-1))
        D[ind[i],ind[i]-q2] <- 1
  D[((q1-1)*q2+1),((q1-1)*q2+1-q2)] <- 1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D  ## t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D  ## t(D)%*%S[[2]]%*%D

  b <- rep(1,q1*q2-1)  # define vector of 0's & 1's for model parameters identification
  b[ind] <- 0
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmd2" specific stuff ...
  object$knots <- list()
  if (is.null(xk))
     object$knots[[1]] <- rep(0,0,0)
  else object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tesmd2.smooth"  # Give object a class
  object
}



####################################################################

marginal.matrices.tesmi2.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 1st unconstrained smooth
{
  # get marginal model matrix for the first unconstrained smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second monotonic smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)

  # create the penalty matrix...
  S <- list()
  # get penalty matrix for the first unconstrained smooth...
  I2 <- diag(q2)
  P <- diff(diag(q1),difference=1)
  S[[1]]<- P %x% I2
  c1 <- c(1,-2,1)
  c2 <- c(1,rep(0,q2-1))
  c <- c1 %x% c2
  for (i in 1:(q1-2)){
       S[[1]][q2*(i-1)+1,(q2*(i-1)+1):(q2*(i-1)+length(c))] <- c
  }
  i <- q1-1
  S[[1]][q2*(i-1)+1,(q2*(i-1)+1):ncol(S[[1]])] <- rep(0,2*q2)

  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd monotonic smooth...
  I2 <- diff(diag(q2-1),difference=1) 
  P <- matrix(0,q2-1,q2)
  P[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%P
  S[[2]] <- crossprod(S[[2]]) ## t(S[[2]])%*%S[[2]]

 list(X1=X1, X2=X2, S=S)
}


###########################################################################
## Prediction matrix for the `tesmd2` smooth class *************************


Predict.matrix.tesmd2.smooth<-function(object,data)
## prediction method function for the `tesmd2' smooth class
{ if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
  IS[1:q2,1]<-1
  for (j in 2:q2)  IS[j,2:j] <- -1
  I <- diag(q1)
  Sig <- I%x%IS
  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}


############################################################################ 
## Tensor product P-spline construction with single monotone increasing   ##
## constraint wrt the first covariate                                     ##
############################################################################


smooth.construct.tesmi1.smooth.spec<- function(object, data, knots)
## construction of the double monotone increasing smooth surface
{ 
 # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      {m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
       object$p.order <- m
  }
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  object$p.order[is.na(object$p.order)] <- 2
  if (object$bs.dim[1]==-1) { # set the default values for q1 and q2
      q1 <- object$bs.dim[1] <- 7
      q2 <- object$bs.dim[2] <- 7
  }
  else if (length(object$bs.dim)==1){
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
      object$bs.dim <- rep(object$bs.dim, 2)
  }
  else {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
    
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
  
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
  {      zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
   }
  if (length(xk)!=nk1 ) # right number of knots?
      stop(paste("there should be ",nk1," supplied knotsfor the x"))
  if (length(zk)!=nk2) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  
  #  get model matrix-------------
  # get marginal model matrices and penalties... 
  bm <- marginal.matrices.tesmi1.ps(x,z,xk,zk,m,q1,q2)
 
  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S      
  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
  IS[1:q1,1]<-1
  for (j in 2:q1)  IS[j,2:j] <- 1
  I <- diag(q2)
  Sig <- IS%x%I

  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D <- diag(q1*q2)
  D <- D[,-q2]
  D1 <- t(diff(diag(q2)))
  D[1:q2,1:(q2-1)] <- D1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D ## t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D ## t(D)%*%S[[2]]%*%D

  b <- rep(1,q1*q2-1)  # define vector of 0's & 1's for model parameters identification
  b[1:(q2-1)] <- rep(0, q2-1) 
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmi1" specific stuff ...
  object$knots <- list()
  object$knots[[1]] <- xk
  if (is.null(zk))
     object$knots[[2]] <- rep(0,0,0)
  else object$knots[[2]] <- zk
  object$m <- m
       
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tesmi1.smooth"  # Give object a class
  object
}



####################################################################

marginal.matrices.tesmi1.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 2nd unconstrained smooth
{
  # get marginal model matrix for the first monotonic smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second unconstrained smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)
   # create the penalty matrix...
  S <- list()
  # get the penalty matrix for the first monotone smooth...
  I2<- diag(q2)
  P <- diff(diag(q1-1),difference=1) 
  Pm1 <- matrix(0,q1-1,q1) # marginal sqrt penalty
  Pm1[2:(q1-1),2:q1] <- P
  S[[1]]<- Pm1%x%I2
  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd smooth
  I2 <- diff(diag(q2),difference=1) 
  I21<- diff(diag(q2),difference=2)
  I1 <- diag(q1)
  S[[2]] <-matrix(0,q2-2+(q1-1)*(q2-1), q1*q2)
  S[[2]][1:(q2-2),] <- t(I1[1,])%x%I21
  S[[2]][(q2-1):nrow(S[[2]]),] <- I1[2:q1,]%x%I2
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]
 list(X1=X1, X2=X2, S=S)
}


## Prediction matrix for the `tesmi1` smooth class *************************


Predict.matrix.tesmi1.smooth<-function(object,data)
## prediction method function for the `tesmi1' smooth class
{ if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
  }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q1,q1)   # Define marginal matrix of Sigma
  IS[1:q1,1]<-1
  for (j in 2:q1)  IS[j,2:j] <- 1
  I <- diag(q2)
  Sig <- IS%x%I

  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}



############################################################################ 
## Tensor product P-spline construction with single monotone decreasing   ##
## constraint wrt the second covariate                                    ##
############################################################################


smooth.construct.tesmi2.smooth.spec<- function(object, data, knots)
## construction of the double monotone increasing smooth surface
{ 
  # require(splines)
  if (object$dim !=2)
      stop("the number of covariates should be two")
  if (length(object$p.order)==1)
      { m <- rep(object$p.order, 2) # if a single number is supplied the same
             ## order of P-splines is provided for both marginal smooths
        object$p.order <- m
  }
  else m <- object$p.order
  m[is.na(m)] <- 2  # the default order is 2 (cubic P-spline)
  object$p.order[is.na(object$p.order)] <- 2
  if (object$bs.dim[1]==-1) { # set the default values for q1 and q2
      q1 <- object$bs.dim[1] <- 7
      q2 <- object$bs.dim[2] <- 7
  }
  else if (length(object$bs.dim)==1){
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
      object$bs.dim <- rep(object$bs.dim, 2)
  }
  else {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  if (is.na(q1)) q1 <- object$bs.dim[1] <- 7  # the default basis dimension is 7
  if (is.na(q2)) q2 <- object$bs.dim[2] <- 7
  
  nk2 <- q2+m[2]+2 ## number of knots for the 2nd smooth
  nk1 <- q1+m[1]+2 ## number of knots for the 1st smooth
  if (nk1<=0 || nk2<=0) stop("either k[1] or k[2] too small for m")
          
  ## the values of the first covariate...   
  x <- data[[object$term[1]]]  
  xk <- knots[[object$term[1]]] ## will be NULL if none supplied
  z <- data[[object$term[2]]]  ## the values of the second covariate
  zk <- knots[[object$term[2]]] ## will be NULL if none supplied
  if (is.null(xk)) # space knots through the values of the 1st covariate
  { n<-length(x)
    xk<-rep(0,q1+m[1]+2)
    xk[(m[1]+2):(q1+1)]<-seq(min(x),max(x),length=q1-m[1])
    for (i in 1:(m[1]+1)) {xk[i]<-xk[m[1]+2]-(m[1]+2-i)*(xk[m[1]+3]-xk[m[1]+2])}
    for (i in (q1+2):(q1+m[1]+2)) {xk[i]<-xk[q1+1]+(i-q1-1)*(xk[m[1]+3]-xk[m[1]+2])}
    knots[[object$term[1]]] <- xk
  }
  if (n != length(z))
     stop ("arguments of smooth not same dimension")
  if (is.null(zk)) # space knots through the values of the 2nd covariate
       { zk<-rep(0,q2+m[2]+2)
         zk[(m[2]+2):(q2+1)]<-seq(min(z),max(z),length=q2-m[2])
         for (i in 1:(m[2]+1)) {zk[i]<-zk[m[2]+2]-(m[2]+2-i)*(zk[m[2]+3]-zk[m[2]+2])}
         for (i in (q2+2):(q2+m[2]+2)) {zk[i]<-zk[q2+1]+(i-q2-1)*(zk[m[2]+3]-zk[m[2]+2])}
         knots[[object$term[2]]] <- zk
  }
  if (length(zk)!=nk2 ) # right number of knots?
      stop(paste("there should be ",nk2," supplied knots for z"))
  if (length(xk)!=nk1) # right number of knots?
      stop(paste("there should be ",nk1," supplied knots for x"))
  
  #  get model matrix-------------
  # get marginal model matrices and penalties... 
  bm <- marginal.matrices.tesmi2.ps(x,z,xk,zk,m,q1,q2)
  
  X1 <- bm$X1
  X2 <- bm$X2
  S <- bm$S  
  # get the overall model matrix...
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- X1[i,]%x%X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
  
  IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
  IS[1:q2,1]<-1
  for (j in 2:q2)  IS[j,2:j] <- 1
  I <- diag(q1)
  Sig <- I%x%IS
    
  # get model matrix
  X <- X%*%Sig
  # apply identifiability constraint 
  D<- diag(q1*q2)
  D<-D[,-((q1-1)*q2+1)]
  ind <- rep(0,q1-1) # get index number for the cells to be changed to "-1"
  for (i in 1:(q1-1)){
        ind[i] <- (i-1)*q2+1
        D[ind[i],ind[i]] <- -1
      }
  for (i in 2:(q1-1))
        D[ind[i],ind[i]-q2] <- 1
  D[((q1-1)*q2+1),((q1-1)*q2+1-q2)] <- 1
  X <- X%*%D 

  object$X <- X # the finished model matrix with identifiability constraint
 
  # create the penalty matrix
  object$S <- list()
  object$S[[1]] <- crossprod(D,S[[1]])%*%D  ## t(D)%*%S[[1]]%*%D
  object$S[[2]] <- crossprod(D,S[[2]])%*%D  ## t(D)%*%S[[2]]%*%D

  b <- rep(1,q1*q2-1)  # define vector of 0's & 1's for model parameters identification
  b[ind] <- 0
  object$p.ident <- b  
  object$rank <- ncol(object$X)-1  # penalty rank
  object$null.space.dim <- m+1  # dim. of unpenalized space
  object$C <- matrix(0, 0, ncol(X)) # to have no other constraints 
  object$Zc <- D   # identifiability constraint matrix

  ## store "tesmi2" specific stuff ...
  object$knots <- list()
  if (is.null(xk))
     object$knots[[1]] <- rep(0,0,0)
  else object$knots[[1]] <- xk
  object$knots[[2]] <- zk
  object$m <- m
   
  object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  class(object)<-"tesmi2.smooth"  # Give object a class
  object
}



####################################################################

marginal.matrices.tesmi2.ps <- function(x,z,xk,zk,m,q1,q2)
## function to get marginal model matrices and penalties in the overall 
## model coefficients in case of P-splines basis for the 1st unconstrained smooth
{
  # get marginal model matrix for the first unconstrained smooth... 
  X1 <- splineDesign(xk,x,ord=m[1]+2)
  # get marginal model matrix for the second monotonic smooth...
  X2 <- splineDesign(zk,z,ord=m[2]+2)

  # create the penalty matrix...
  S <- list()
  # get penalty matrix for the first unconstrained smooth...
  I2 <- diag(q2)
  P <- diff(diag(q1),difference=1)
  S[[1]]<- P %x% I2
  c1 <- c(1,-2,1)
  c2 <- c(1,rep(0,q2-1))
  c <- c1 %x% c2
  for (i in 1:(q1-2)){
       S[[1]][q2*(i-1)+1,(q2*(i-1)+1):(q2*(i-1)+length(c))] <- c
  }
  i <- q1-1
  S[[1]][q2*(i-1)+1,(q2*(i-1)+1):ncol(S[[1]])] <- rep(0,2*q2)

  S[[1]] <- crossprod(S[[1]])  ## t(S[[1]])%*%S[[1]]

  # get penalty for the 2nd monotonic smooth...
  I2 <- diff(diag(q2-1),difference=1) 
  P <- matrix(0,q2-1,q2)
  P[2:(q2-1),2:q2] <- I2  # marginal sqrt penalty
  I1 <- diag(q1)
  S[[2]] <- I1%x%P
  S[[2]] <- crossprod(S[[2]])  ## t(S[[2]])%*%S[[2]]

 list(X1=X1, X2=X2, S=S)
}



############################################################
## Prediction matrix for the `tesmi2` smooth class ....   ##
############################################################


Predict.matrix.tesmi2.smooth<-function(object,data)
## prediction method function for the `tesmi2' smooth class
{ if (length(object$bs.dim)==1)
      q1 <- q2 <- object$bs.dim # if `k' is supplied as a single number, the same
             ## basis dimension is provided for both marginal smooths
  else  {q1 <- object$bs.dim[1]; q2 <- object$bs.dim[2]}
  
  bm <- marginal.linear.extrapolation(object, data)
  n <- length(data[[object$term[1]]])
  X <- matrix(0,n,q1*q2)  # model matrix
  for (i in 1:n)
    {  X[i,] <- bm$X1[i,] %x% bm$X2[i,] # Kronecker product of two rows of marginal model matrices
    }
  # get a matrix Sigma -----------------------
  IS <- matrix(0,q2,q2)   # Define submatrix of Sigma
  IS[1:q2,1]<-1
  for (j in 2:q2)  IS[j,2:j] <- 1
  I <- diag(q1)
  Sig <- I%x%IS
    
  # get final model matrix
  X <- X%*%Sig
  X  # return the prediction matrix
}


