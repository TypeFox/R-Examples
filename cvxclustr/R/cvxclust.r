#' Convex Clustering Path via Variable Splitting Methods
#' 
#' \code{cvxclust} estimates the convex clustering path via variable splitting methods: ADMM and AMA. This function
#' is a wrapper function that calls either \code{cvxclust_path_admm} or \code{cvxclust_path_ama} (the default) to perform the computation.
#' Required inputs include a data matrix \code{X} (rows are features; columns are samples), a vector of weights
#' \code{w}, and a sequence of regularization parameters \code{gamma}.
#' Two penalty norms are currently supported: 1-norm and 2-norm. Both ADMM and AMA admit acceleration schemes at little additional computation.
#' Acceleration is turned on by default.
#' 
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param method Algorithm to use: "admm" or "ama"
#' @param gamma A sequence of regularization parameters.
#' @param nu A positive penalty parameter for quadratic deviation term.
#' @param tol The convergence tolerance.
#' @param max_iter The maximum number of iterations.
#' @param type An integer indicating the norm used: 1 = 1-norm, 2 = 2-norm.
#' @param accelerate If \code{TRUE} (the default), acceleration is turned on.
#' @return \code{U} A list of centroid matrices.
#' @return \code{V} A list of centroid difference matrices.
#' @return \code{Lambda} A list of Lagrange multiplier matrices.
#' @export
#' @author Eric C. Chi, Kenneth Lange
#' @seealso \code{\link[cvxclustr]{cvxclust_path_ama}} and \code{\link[cvxclustr]{cvxclust_path_admm}} for estimating the clustering path with AMA or ADMM. 
#' \code{\link[cvxclustr]{kernel_weights}} and \code{\link[cvxclustr]{knn_weights}} compute useful weights. 
#' To extract cluster assignments from the clustering path use \code{\link[cvxclustr]{create_adjacency}} and \code{\link[cvxclustr]{find_clusters}}.
#' @examples
#' ## Clusterpaths for Mammal Dentition
#' data(mammals)
#' X <- as.matrix(mammals[,-1])
#' X <- t(scale(X,center=TRUE,scale=FALSE))
#' n <- ncol(X)
#' 
#' ## Pick some weights and a sequence of regularization parameters.
#' k <- 5
#' phi <- 0.5
#' w <- kernel_weights(X,phi)
#' w <- knn_weights(w,k,n)
#' gamma <- seq(0.0,43, length.out=100)
#' 
#' ## Perform clustering
#' sol <- cvxclust(X,w,gamma)
#' 
#' ## Plot the cluster path
#' library(ggplot2)
#' svdX <- svd(X)
#' pc <- svdX$u[,1:2,drop=FALSE]
#' pc.df <- as.data.frame(t(pc)%*%X)
#' nGamma <- sol$nGamma
#' df.paths <- data.frame(x=c(),y=c(), group=c())
#' for (j in 1:nGamma) {
#'   pcs <- t(pc)%*%sol$U[[j]]
#'   x <- pcs[1,]
#'   y <- pcs[2,]
#'   df <- data.frame(x=pcs[1,], y=pcs[2,], group=1:n)
#'   df.paths <- rbind(df.paths,df)
#' }
#' X_data <- as.data.frame(t(X)%*%pc)
#' colnames(X_data) <- c("x","y")
#' X_data$Name <- mammals[,1]
#' data_plot <- ggplot(data=df.paths,aes(x=x,y=y))
#' data_plot <- data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
#' data_plot <- data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name),
#'   position=position_jitter(h=0.125,w=0.125))
#' data_plot <- data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
#' data_plot <- data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
#' data_plot + theme_bw()
#' 
#' ## Output Cluster Assignment at 10th gamma
#' A <- create_adjacency(sol$V[[10]],w,n)
#' find_clusters(A)
#' 
#' ## Visualize Cluster Assignment
#' G <- graph.adjacency(A, mode = 'upper')
#' plot(G,vertex.label=as.character(mammals[,1]),vertex.label.cex=0.65,vertex.label.font=2)
cvxclust <- function(X,w,gamma,method="ama",nu=1,tol=1e-3,max_iter=1e4,type=2,accelerate=TRUE) {
  if (!is.null(method) && !(method %in% c("ama","admm")))
    stop("method must be 'ama', 'admm', or NULL.")
  if (method == "ama") {
    cvxclust_obj <- cvxclust_path_ama(X,w,gamma,nu=nu,tol=tol,max_iter=max_iter,type=type,accelerate=accelerate)
  } else {
    cvxclust_obj <- cvxclust_path_admm(X,w,gamma,nu=nu,tol_abs=tol,tol_rel=tol,max_iter=max_iter,type=type,accelerate=accelerate)
  }

  return(cvxclust_obj)
}
