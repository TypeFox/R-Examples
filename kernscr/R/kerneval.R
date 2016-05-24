#'Evaluation of kernels
#'
#'@param tZ a \code{P x N} matrix of genomic covariates (i.e., the usual data array Z transposed)
#'
#'@param kernel which kernel is evaluated by \code{kerneval}. Possible values include currently implemented kernels
#'designated by a character string \code{"linear"}, \code{"poly"} and \code{"gaussian"}.
#'Otherwise can also be a user-defined function (see \code{kernel_func}).
#'
#'@param ... other arguments to be passed to be passed to the evaluated kernel function.
#'
#'@details \code{kernelEval} works only for gaussian, polynomial and linear kernels currently.
#'
#'@return \code{kernelEval}, \code{linKernelEval}, \code{gaussKernelEval}, and \code{genericKernelEval}
#'return an \code{N x N} matrix with entries \code{K(Z[i,], Z[j,])} [persons i,j]
#'
#'@rdname kernelEval
#'@keywords internal
#'@export
kernelEval <- function(tZ, kernel = c("linear", "poly", "gaussian"), ... ){
  if (kernel  == "linear"){
    return(linKernelEval(tZ))
  }
  else if (kernel == "poly"){
    return(polyKernelEval(tZ, ...))
  }
  else if(kernel == "gaussian"){
    return(gaussKernelEval(tZ, ...))
  }
  else{
    return(genericKernelEval(tZ, kernel, ...))
  }
}


#'@rdname kernelEval
#'@keywords internal
linKernelEval <- function(tZ){
  tZ <- as.matrix(tZ) # this is the size Gram matrix
  return(t(tZ)%*%tZ)
}


#'@param sigma standard-deviation parameter for the \code{"gaussian"} kernel.
#'@rdname kernelEval
#'@keywords internal
gaussKernelEval <- function(tZ, sigma = 1){

  tZ <- as.matrix(tZ) # this is the size Gram matrix
  m <- ncol(tZ)

  # computing the Gram matrix
  Z_rep <- tZ[,rep(1:m, rep(m, m))] - tZ[,rep(1:m, m)]
  return(matrix(exp(-colSums(Z_rep * Z_rep)/(2*sigma^2)), m, byrow = T))
}

#'@rdname kernelEval
#'@param a TODO of the polynomial for the \code{"poly"}. Default is \code{0}
#'@param d power of the polynomial.  Default is \code{2} (quadratic kernel).
#'@keywords internal
polyKernelEval <- function(tZ, a = 0, d = 2){
  tZ <- as.matrix(tZ) # this is the size Gram matrix
  m <- ncol(tZ)
  return((matrix(as.vector(t(tZ)%*%tZ) + rep(a,rep(m*m)), ncol = m*m, byrow = T))^d)
}




#'@rdname kernelEval
#'@param kernel_func a function, whose first argument should be \code{tZ}
#'@details \code{genericKernelEval}
genericKernelEval <- function(tZ, kernel_func, ...){
  ## this function is too slow; don't know how to optimize it as of today...
  ## maybe somekind of apply to the vector columns
  if (!is.function(kernel_func)){
    cat("kernel_func must be a kernel function defined for the columns of X")
    return(0)
  }

  Z <- as.matrix(tZ) # this is the size Gram matrix
  m <- ncol(tZ)

  # computing the Gram matrix
  G <- matrix(data = 0, nrow = m, ncol = m)
  for (i in 1:m){
    for (j in 1:m){
      G[i,j] <- kernel_func(as.vector(tZ[,i]), as.vector(tZ[,j]), ...)
    }
  }

  return(G)
}









#'@rdname kernelEval
#'@param rho either a single rho to evaluate the kernel at, or a vector of rhos
#'@return\code{gaussKernelEval_multipleRhos} and \code{polyKernelEval_multipleRhos} return
#'a matrix of dimension \code{Q x N^2}, where \code{Q} is the \code{length} of \code{rho},
#'each row corresponds to a rho (puns!) to get the actual kernel matrix associated with a particular
#'value of rho, if output is \code{G}, take \code{matrix(G[i,], N)}
#'@keywords internal
#'@export
gaussKernelEval_multipleRhos <- function(tZ, rho){
  tZ <- as.matrix(tZ) # this is the size Gram matrix
  m <- ncol(tZ)

  # computing the Gram matrix
  tZ_rep <- tZ[,rep(1:m, rep(m, m))] - tZ[,rep(1:m, m)]
  return(exp(-matrix(colSums(tZ_rep*tZ_rep)[rep(1:ncol(tZ_rep), length(rho))], nrow=length(rho), byrow=T)/(rho)))
}

#'@rdname kernelEval
#'@keywords internal#'
#'@details For \code{polyKernelEval_multipleRhos}, one should have \code{rho} > 0 to get
#'basis of monomials up to degree \code{d}
#'
#'@export

polyKernelEval_multipleRhos <- function(tZ, rho, d=2){
  tZ <- as.matrix(tZ) # this is the size Gram matrix
  m <- ncol(tZ)
  return((matrix(as.vector(t(tZ)%*%tZ) + rep(rho, rep(m*m,length(rho))), ncol = m*m, byrow = TRUE))^d)
}



