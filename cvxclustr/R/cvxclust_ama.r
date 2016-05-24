#' Convex clustering via AMA
#' 
#' \code{cvxclust_ama} performs convex clustering via AMA. This is an R wrapper function around C code.
#' Dimensions of various arguments are as follows:
#' \itemize{
#' \item{n is the number of data points}
#' \item{p is the number of features}
#' \item{k is the number non-zero weights.}
#' }
#' Note that the indices matrices 'M1', 'M2', and 'ix' take on values starting at 0 to match the indexing conventions of C.
#' 
#' @param X The p-by-n data matrix whose columns are to be clustered.
#' @param Lambda The p-by-k matrix of Lagrange multipliers.
#' @param ix The k-by-2 matrix of index pairs.
#' @param M1 Index set used to track nonzero weights.
#' @param M2 Index set used to track nonzero weights.
#' @param s1 Index set used to track nonzero weights.
#' @param s2 Index set used to track nonzero weights.
#' @param w A vector of k positive weights.
#' @param gamma The regularization parameter controlling the amount of shrinkage.
#' @param nu The initial step size parameter when backtracking is applied. Otherwise it is a fixed step size in which case there are no guarantees of convergence if it exceeds \code{2/ncol(X)}.
#' @param type An integer indicating the norm used: 1 = 1-norm, 2 = 2-norm.
#' @param max_iter The maximum number of iterations.
#' @param tol The convergence tolerance.
#' @param accelerate If \code{TRUE} (the default), acceleration is turned on.
#' @return \code{U} A list of centroid matrices.
#' @return \code{V} A list of centroid difference matrices.
#' @return \code{Lambda} A list of Lagrange multiplier matrices.
#' @return \code{nu} The final step size used.
#' @return \code{primal} The primal objective evaluated at the final iterate.
#' @return \code{dual} The dual objective evaluated at the final iterate.
#' @return \code{iter} The number of iterations taken.
#' @export
#' @author Eric C. Chi, Kenneth Lange
#' @useDynLib cvxclustr
#' @examples
#' ## Create random problem
#' seed <- 12345
#' p <- 10
#' n <- 20
#' rnd_problem <- create_clustering_problem(p,n,seed)
#' X <- rnd_problem$X
#' ix <- rnd_problem$ix
#' M1 <- rnd_problem$M1
#' M2 <- rnd_problem$M2
#' s1 <- rnd_problem$s1
#' s2 <- rnd_problem$s2
#' w  <- rnd_problem$w
#' nK <- nrow(ix)
#' Lambda <- matrix(rnorm(p*nK),p,nK)    
#' gamma <- 0.1
#' nu <- 1.999/n
#' max_iter <- 1e6
#' tol <- 1e-15
#' sol_ama <- cvxclust_ama(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,max_iter=max_iter,tol=tol)
cvxclust_ama <- function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,type=2,max_iter=1e2,tol=1e-4,accelerate=TRUE) {
  p <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  if (!is.null(type) && !(type %in% c(1,2)))
    stop("type must be 1, 2, or NULL. Only 1-norm and 2-norm penalties are currently supported.")
  if (nu >= 2/n) {
    warning("The stepsize nu may be too large. Setting it to maximum of 1/n and 1/AnM.")
    nu <- max(AMA_step_size(w,n),1/n)
  }
  nK <- as.integer(ncol(Lambda))
  mix1 <- as.integer(nrow(M1))
  mix2 <- as.integer(nrow(M2))
  storage.mode(X) <- "double"
  storage.mode(Lambda) <- "double"
  U <- matrix(0,p,n)
  storage.mode(U) <- "double"
  V <- matrix(0,p,nK)
  storage.mode(V) <- "double"
  storage.mode(ix) <- "integer"
  storage.mode(M1) <- "integer"
  storage.mode(M2) <- "integer"
  s1 <- as.integer(s1)
  s2 <- as.integer(s2)
  w <- as.double(w)
  gamma <- as.double(gamma)
  nu <- as.double(nu)
  max_iter <- as.integer(max_iter)
  tol <- as.double(tol)
  primal <- double(max_iter)
  dual <- double(max_iter)
  if (accelerate) {
    fxname <- 'convex_cluster_ama_acc'
  } else {
    fxname <- 'convex_cluster_ama'
  }  
  
  type <- as.integer(type)
  sol <- .C(fxname,X=X,Lambda=Lambda,U=U,V=V,p=p,n=n,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,type=type,
                 s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}
