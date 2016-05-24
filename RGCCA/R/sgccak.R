#' The function sgccak() is called by sgcca() and does not have to be used by the user. 
#' sgccak() enables the computation of SGCCA block components, outer weight vectors, etc., 
#' for each block and each dimension. 
#' @param A  A list that contains the \eqn{J} blocks of variables from which block components are constructed. 
#' It could be eiher the original matrices (\eqn{\mathbf{X}_1, \mathbf{X}_2, \ldots, \mathbf{X}_J}) or the residual matrices (\eqn{\mathbf{X}_{h1}, \mathbf{X}_{h2}, \ldots, \mathbf{X}_{hJ}}).
#' @param C  A design matrix that describes the relationships between blocks.
#' @param c1 A \eqn{1 \times J} vector that contains the value of c1 applied to each block. The L1 bound on a[[j]] is 
#' \deqn{ \|a_{j}\|_{\ell_1} \leq c_1[j] \sqrt{p_j}.}
#' with \eqn{p_j} the number of variables of \eqn{\mathbf{X}_j} and with c1[j] between 0 and 1 (larger L1 bound corresponds to less penalization).
#' @param scheme  Either "horst", "factorial" or "centroid" (default: centroid).
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param init Mode of initialization of the SGCCA algorithm. Either by Singular Value Decompostion ("svd") or random ("random") (default: "svd").
#' @param bias Logical value for biaised (\eqn{1/n}) or unbiaised (\eqn{1/(n-1)}) estimator of the var/cov.
#' @param verbose  Reports progress while computing, if verbose = TRUE (default: TRUE).
#' @param tol Stopping value for convergence.
#' @return \item{Y}{A \eqn{n \times J} matrix of SGCCA block components.}
#' @return \item{a}{A list of \eqn{J} elements. Each element contains the outer weight vector of each block.}
#' @return \item{crit}{The values of the objective function at each iteration of the iterative procedure.}
#' @return \item{converg}{Speed of convergence of the alogrithm to reach the tolerance.}
#' @return \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @return \item{C}{A design matrix that describes the relationships between blocks (user specified).}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @title Internal function for computing the SGCCA parameters (SGCCA block components, outer weight vectors etc.)
#' @export sgccak
sgccak <-  function(A, C, c1 = rep(1, length(A)), scheme = "centroid", scale = FALSE,
                    tol = .Machine$double.eps, init="svd", bias = TRUE, verbose = TRUE){

  J <- length(A)
  pjs = sapply(A, NCOL)
  AVE_X <- rep(0, J)
  # Data standardization 
  if (scale == TRUE) A <- lapply(A, function(x) scale2(x, bias = bias))
  #  Choose J arbitrary vectors
  if (init=="svd") {
    #SVD Initialisation of a_j or \alpha_j
    a <- lapply(A, function(x) return(svd(x,nu=0,nv=1)$v)) #
  } else if (init=="random") {
    a <- lapply(pjs,rnorm)
  } else {
    stop("init should be either random or svd.")
  }
  if (any( c1 < 0 | c1 > 1 )) stop("L1 constraints must be between 0 and 1.") 
    
  const <- c1*sqrt(pjs)
  #	Apply the constraints of the general otpimization problem
  #	and compute the outer components
  iter <- 1
  converg <- crit <- numeric()
  Y <- Z <- matrix(0,NROW(A[[1]]),J)
  for (q in 1:J){
      Y[,q] <- apply(A[[q]],1,miscrossprod,a[[q]])
      a[[q]] <- soft.threshold(a[[q]], const[q])   
      a[[q]] <- as.vector(a[[q]])/norm2(a[[q]])
  }	 
  a_old <- a
  g <- function(x) switch(scheme,horst=x,factorial=x**2,centroid=abs(x))
  
  
  repeat{
    Yold <- Y
    for (q in 1:J){
      if (scheme == "horst"    ) CbyCovq <- C[q,]
      if (scheme == "factorial") CbyCovq <- C[q,]*cov2(Y, Y[, q], bias = bias)
      if (scheme == "centroid" ) CbyCovq <- C[q,]*sign(cov2(Y, Y[,q], bias = bias))
      Z[,q] <- rowSums(mapply("*",CbyCovq,as.data.frame(Y)))
      a[[q]] <- apply( t(A[[q]]),1,miscrossprod, Z[,q])
      a[[q]] <- soft.threshold(a[[q]], const[q])
      a[[q]] <- as.vector(a[[q]])/norm2(a[[q]]) 
      Y[,q] <- apply(A[[q]], 1, miscrossprod,a[[q]])
    }
    
    num_converg <- sum((rowSums(Yold) - rowSums(Y))^2)
    den_converg <- sum(rowSums(Yold)^2)
    converg[iter] <- num_converg/den_converg
    #check for convergence of the SGCCA alogrithm to a fixed point of the stationnary equations
    stationnary_point <- rep(FALSE, J)
    for (ind in 1:J) 
      stationnary_point[ind] <- sum(round(abs(a_old[[ind]]-a[[ind]]), 8) < tol) == NCOL(A[[ind]])
    a_old <- a
    crit[iter] <- sum(C*g(cov2(Y, bias = bias)))
    if (iter > 1000) warning("The SGCCA algorithm did not converge after 1000 iterations.")
    if ((converg[iter] < tol & sum(stationnary_point) == J) | iter > 1000)  break
    iter <- iter + 1 
  }
  if(sum(stationnary_point) == J & verbose) cat("The SGCCA algorithm converged to a fixed point of the stationary equations after", iter-1, "iterations \n") 
  if (verbose) plot(crit, xlab = "iteration", ylab = "criteria")
  
  for (q in 1:J) if(sum(a[[q]]!=0) <= 1) warning(sprintf("Deflation failed because only one variable was selected for block #",q))
  
  AVE_inner  <- sum(C*cor(Y)^2/2)/(sum(C)/2) # AVE inner model
            
  result <- list(Y = Y, a = a, crit = crit[which(crit != 0)], 
                converg = converg[which(converg != 0)], 
                AVE_inner = AVE_inner, C = C, c1, scheme = scheme)
  return(result)                                                                          
}


