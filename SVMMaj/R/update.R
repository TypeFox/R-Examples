getUpdate <- function( X,kernel =  vanilladot())
{
X <-data.matrix(X)
#===================================================================
#DETERMINE WHICH METHOD TO USE
#-------------------------------------------------------------------
#LINEAR KERNEL: QR-DECOMPOSITION
# [ X = Z Q' ]
if(class(kernel)=='vanillakernel'){
    RQ   <- qr(t(X))
    
    method <- list(type  = 'QR',
                      x  = RQ$rank,
                      X  = X,
                      Z  = t(qr.R(RQ))[order(RQ$pivot),1:RQ$rank],
                      RQ = RQ)   
                             
      
} else {
#NONLINEAR KERNEL: PIVOTED CHOLESKY DECOMPOSITION
# [ K   = Phi Phi' = Z Z' = P R R' P' ]
# [ Phi = P R Q', Q' = UNDEFINED  ]
  options(warn=-1)
    Rp   <- chol(kernelMatrix(kernel,X),pivot=TRUE)
  options(warn=0)
    r    <- sum(diag(Rp)/Rp[1]>1e-3)
    pivot<- order(attr(Rp,'pivot'))
    
    method <- list(type   = 'Cholesky',
                      x   = r,
                      X   = X,
                      Z   = t(Rp[1:r ,order(attr(Rp,'pivot'))]))
}
#gc()
class(method) <- 'theta'
method$kernel <- kernel
method$linear <- class(kernel)=='vanillakernel'
return(method)
}


#===================================================================
#CALCULATE THE INTERCEPT alpha GIVEN THE model AND theta
#-------------------------------------------------------------------
alpha.theta <- function(object,theta) return(theta[1,])

#===================================================================
#CALCULATE THE WEIGHTS beta GIVEN THE model AND theta
#-------------------------------------------------------------------
beta.theta <- function(object,theta) {
    if(!object$linear) return(NULL)
    b <- qr.Q(object$RQ)[,1:object$x] %*% theta[-1,]
    row.names(b)<- colnames(object$X)
    return(b)
}


#===================================================================
#CALCULATE THE PREDICTED VALUE q GIVEN THE model, theta AND x
#-------------------------------------------------------------------
qhat.theta <- function(object,theta,Xnew) {
    if(!object$linear){
        QR   =  qr(object$Z)
        qhat <- theta[1,] +
                kernelMatrix(object$kernel,data.matrix(Xnew),object$X) %*%
                qr.Q(QR) %*% solve(t(qr.R(QR)), theta[-1,])
    } else
        qhat <- theta[1] + data.matrix(Xnew) %*% beta.theta(object,theta)
    return(drop(qhat))
}

