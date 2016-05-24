

#' Hinge Loss function for SVM
#'
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x in case of misprediction. 
#'        Vector length should match length(y), but values are cycled if not of identical size. 
#'        Default to 1 so we define a standard 0/1 loss for SVM classifier. 
#'        The parameter might be useful to adapt SVM learning in case of unbalanced class distribution.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
hingeLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y) || !all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * pmax(0,1-y*f)
    grad <- loss.weights * (loss>0) * (-y)
    val <- sum(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to maximize area under the ROC curve
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
rocLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y) || !all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1')    
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- rep(w,length.out=ncol(x))
    c <- x %*% w - 0.5*y
    o <- order(c)
    
    sp <- cumsum(y[o]==+1)
    sm <- sum(y==-1) - cumsum(y[o]==-1)
    l <- numeric(length(o))
    l[o] <- ifelse(y[o]==-1,sp,-sm)
    l <- l/(sum(y==-1)*sum(y==+1))
    
    val <- crossprod(l,c)
    gradient(val) <- crossprod(l,x)
    return(val)
  }
}



#' The loss function for ordinal regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y integer vector of positive values (>=1) representing the training labels for each instance in x
#' @param C the cost matrix to use, C[i,j] being the cost for predicting label i instead of label j.
#' @param impl either the string "loglin" or "quadratic", that define the implementation to use for the computation of the loss.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#' # -- Load the data
#' x <- data.matrix(iris[1:4])
#' y <- as.integer(iris$Species)
#' 
#' # -- Train the model
#' m <- bmrm(ordinalRegressionLoss(x,y),LAMBDA=0.001,EPSILON_TOL=0.0001)
#' m2 <- bmrm(ordinalRegressionLoss(x,y,impl="quadratic"),LAMBDA=0.001,EPSILON_TOL=0.0001)
#' 
#' # -- plot predictions
#' f <- x %*% m$w
#' f2 <- x %*% m2$w
#' layout(1:2)
#' plot(y,f)
#' plot(f,f2,main="compare predictions of quadratic and loglin implementations")
#' 
#' # -- Compute accuracy
#' ij <- expand.grid(i=seq(nrow(x)),j=seq(nrow(x)))
#' n <- tapply(f[ij$i] - f[ij$j]>0,list(y[ij$i],y[ij$j]),sum)
#' N <- table(y[ij$i],y[ij$j])
#' print(n/N)
ordinalRegressionLoss <- function(x,y,C="0/1",impl=c("loglin","quadratic")) {
  impl <- match.arg(impl)

  # check parameters at first call
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')  
  C <- costMatrix(y,C)
  y <- as.integer(y)
  m <- length(y)
  mi <- tabulate(y,nbins=ncol(C))
  M <- (m*m - sum(mi*mi))/2
  C <- C / M
  
  .loglin <- function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    c <- c(f-0.5,f+0.5)
    o <- order(c)
    
    j <- ((o-1)%%m)+1
    
    l <- matrix(0,2*m,length(mi))
    l[cbind(which(o<=m),y[j[o<=m]])] <- 1
    l <- apply(l,2,cumsum)
    
    u <- matrix(0,2*m,length(mi))
    u[cbind(which(o>m),y[j[o>m]])] <- 1
    u <- mi[col(u)] - apply(u,2,cumsum)
    
    Gu <- t(C)[y[j],] * u
    Gu[col(Gu)>=y[j]] <- 0
    Gl <- C[y[j],] * l
    Gl[col(Gl)<=y[j]] <- 0
    
    v <- ifelse(o<=m,-rowSums(Gu), rowSums(Gl))
    r <- sum(v*c[o])
    g <- matrix(NA,m,2)
    g[cbind(j,1 + (o-1)%/%m)] <- v
    g <- rowSums(g)
    
    gradient(r) <- crossprod(g,x)
    return(r)
  }
  
  .quadratic <- function(w) {
      w <- rep(w,length.out=ncol(x)) 
      f <- x %*% w
      
      # alternative computation in quadratic time for debugging purpose only
      z <- expand.grid(i=factor(1:m),j=factor(1:m))
      z <- z[y[z$i] < y[z$j],]
      z <- z[1+f[z$i]-f[z$j]>0,]
      R <- sum(C[cbind(y[z$i],y[z$j])] * (1+f[z$i]-f[z$j]))
      gradient(R) <- colSums(C[cbind(y[z$i],y[z$j])] * (x[z$i,]-x[z$j,]))
      return(R)
  }
  
  switch(impl,loglin=.loglin,quadratic=.quadratic)
}




#' F beta score loss function
#'
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param beta a numeric value setting the beta parameter is the f-beta score
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
fbetaLoss <- function(x,y,beta=1) {
  
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y) || !all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  
  .fbeta <- function(TP,TN,P,N,beta) {
    beta2 <- beta*beta
    (1+beta2)*TP / (TP+N-TN+beta2*P)
  }
  
  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    o <- order(f,decreasing=TRUE)
    op <- o[y[o]==1]
    on <- rev(o[y[o]==-1])
    
    p <- 2*(sum(f[op]) - cumsum(c(0,f[op])))
    n <- 2*(sum(f[on]) - cumsum(c(0,f[on])))
    R <- outer(seq_along(p),seq_along(n),function(i,j) {
      1 - .fbeta(i-1,j-1,length(op),length(on),beta) - p[i] + n[j]
    })
    
    ij <- arrayInd(which.max(R),dim(R))
    Y <- -y
    Y[op[seq_len(ij[1,1]-1L)]] <- 1L
    Y[on[seq_len(ij[1,2]-1L)]] <- -1L
    
    val <- R[ij]
    gradient(val) <- crossprod(x,Y-y)
    return(val)
  }
}


