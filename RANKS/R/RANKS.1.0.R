#################################
## Graph kernels
#################################

########################################################
#' Method to compute the random walk kernel (Smola and Kondor, 2003)
#' It computes RW = (a-1)I + D^(-1/2) * W * D^(-1/2)
#' where I is the identity matrix, W is the weight matrix of an undirected graph,
#' and  D is a diagonal matrix with D_ii = sum_j W_ij
#' @param W : a square symmetric matrix with positive values
#' @param a : numeric. It is correlated to the probability of remaining at the same vertex. Larger a, larger the probability (def. 2)
#' @return  the one-step random walk kernel matrix
#' @export
#' @docType methods
#' @rdname rw.kernel-methods
setGeneric("rw.kernel", 
                 function(W, a=2) standardGeneric("rw.kernel"));  
		 
#' Method to compute the random walk kernel (Smola and Kondor, 2003)
#' It computes RW = (a-1)I + D^(-1/2) * W * D^(-1/2)
#' where I is the identity matrix, W is the weight matrix of an undirected graph, 
#' and  D is a diagonal matrix with D_ii = sum_j W_ij
#' @param W : a square symmetric matrix with positive values
#' @param a : numeric. It is correlated to the probability of remaining at the same vertex. Larger a, larger the probability (def. 2)
#' @return  the one-step random walk kernel matrix
setMethod("rw.kernel", signature(W="matrix"),
  function(W, a=2) {
   n <- nrow(W);
   if (n != ncol(W))
     stop("first arg must be a square matrix");
   if (a < 1)
     stop("coefficient a must be larger or equal than 1");  
   # dyn.load("ranks.so");
   names.var <- rownames(W);
   # computing (a-1)I
   I <- diag(x=a-1, nrow=n);
   # computing D^(-1/2) * W * D^(-1/2)
   diag.D <- apply(W,1,sum);
   diag.D[diag.D==0] <- Inf;
   inv.sqrt.diag.D <- 1/sqrt(diag.D);
   W <- .C("norm_lapl_graph", as.double(W), as.double(inv.sqrt.diag.D), as.integer(n), PACKAGE = "RANKS")[[1]];
   W <- I + W;
   rownames(W) <- colnames(W) <- names.var;
   return(W);
})


#' Method to compute the random walk kernel (Smola and Kondor, 2003)
#' @param W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param a : numeric. It is correlated to the probability of remaining at the same vertex. Larger a, larger the probability (def. 2)
#' @return The one-step random walk kernel matrix
setMethod("rw.kernel", signature(W="graph"),
function(W, a=2) {
  W <- as(W, "matrix");
  return(rw.kernel(W, a));
})


##################################################################
#' Method to compute a p-step random walk kernel.
#' It computes a p-step random walk kernel RW^p
#' @param RW : a square symmetric matrix representing a one-step random walk kernel
#' @param p : integer. Number of steps (def: p=2)
#' @return The p-step random walk kernel matrix
#' @export
#' @docType methods
#' @rdname p.step.rw.kernel-methods
setGeneric("p.step.rw.kernel", 
                 function(RW, p=2) standardGeneric("p.step.rw.kernel"));

#' Method to compute a p-step random walk kernel.
#' It computes a p-step random walk kernel RW^p
#' @param RW : a square symmetric matrix representing a one-step random walk kernel
#' @param p : integer. Number of steps (def: p=2)
#' @return The p-step random walk kernel matrix
setMethod("p.step.rw.kernel", signature(RW="matrix"),
  function(RW, p=2) {
	if (nrow(RW) != ncol(RW))
      stop("first arg must be a square matrix");
	if (p < 2)
      stop("the number of steps a must be larger or equal than 2");
	Res <- RW;
	for (i in 2:p)
      Res <- Res %*% RW;
	return(Res);
})

#' Method to compute a p-step random walk kernel.
#' It computes a p-step random walk kernel RW^p
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param p : integer. Number of steps (def: p=2)
#' @return The p-step random walk kernel matrix
setMethod("p.step.rw.kernel", signature(RW="graph"),
  function(RW, p=2) {
	W <- rw.kernel(RW);
	return(p.step.rw.kernel(W,p));
})

#' Function to compute a identity  kernel from a square symmetric matrix with positive values
#' @param W : a square symmetric matrix with positive values 
#' @param a : unused parameter, maintained for compatibility reasons 
#' @return The same matrix W
identity.kernel <- function(W, a=1) {
   return(W);
}
		 
#' Function to compute a linear  kernel from a feature matrix 
#' It computes K = W  * t(W), where W is a feature matrix
#' @param W : a feature matrix, Rows are examples and columns are features    
#' @param  a : unused parameter, maintained for compatibility reasons     
#' @return A linear kernel matrix 
linear.kernel <- function(W, a=1) {
   return(W %*% t(W));
}

#' Function to compute a gaussian  kernel from a feature matrix 
#' @param W : a feature matrix, Rows are examples and columns are features      
#' @param sigma : a real value representing the sigma parameter of the gaussian  
#' @return A gaussian kernel matrix 
gaussian.kernel <- function(W, sigma=1) {
   n <- nrow(W);
   m <- ncol(W);
   nn <- rownames(W);
   K <- matrix(numeric(n*n), nrow=n);
   K <- .C("gaussian_kernel", as.double(K), as.double(W), as.double(sigma), as.integer(n), as.integer(m), PACKAGE = "RANKS")[[1]];
   K <- matrix(K, nrow=n);
   rownames(K) <- colnames(K) <- nn;
   return(K);
}

#' Function to compute a laplacian  kernel from a feature matrix 
#' @param W : a feature matrix, Rows are examples and columns are features      
#' @param sigma : a real value representing the sigma parameter of the laplacian  
#' @return A laplacian kernel matrix 
laplacian.kernel <- function(W, sigma=1) {
   n <- nrow(W);
   m <- ncol(W);
   nn <- rownames(W);
   K <- matrix(numeric(n*n), nrow=n);
   K <- .C("laplacian_kernel", as.double(K), as.double(W), as.double(sigma), as.integer(n), as.integer(m), PACKAGE = "RANKS")[[1]];
   K <- matrix(K, nrow=n);
   rownames(K) <- colnames(K) <- nn;
   return(K);
}

#' Function to compute a Cauchy  kernel from a feature matrix 
#' @param W : a feature matrix, Rows are examples and columns are features      
#' @param sigma : a real value representing the sigma parameter of the cauchy kernel 
#' @return A Cauchy  kernel matrix 
cauchy.kernel <- function(W, sigma=1) {
   n <- nrow(W);
   m <- ncol(W);
   nn <- rownames(W);
   K <- matrix(numeric(n*n), nrow=n);
   K <- .C("cauchy_kernel", as.double(K), as.double(W), as.double(sigma), as.integer(n), as.integer(m), PACKAGE = "RANKS")[[1]];
   K <- matrix(K, nrow=n);
   rownames(K) <- colnames(K) <- nn;
   return(K);
}

#' Function to compute an inverse multiquadric  kernel from a feature matrix 
#' @param W : a feature matrix, Rows are examples and columns are features      
#' @param v : constant factor (def. 1). v must be larger than 0.
#' @return An inverse multiquadric  kernel matrix 
inv.multiquadric.kernel <- function(W, v=1) {
   n <- nrow(W);
   m <- ncol(W);
   nn <- rownames(W);
   K <- matrix(numeric(n*n), nrow=n);
   K <- .C("inv_mq_kernel", as.double(K), as.double(W), as.double(v), as.integer(n), as.integer(m), PACKAGE = "RANKS")[[1]];
   K <- matrix(K, nrow=n);
   rownames(K) <- colnames(K) <- nn;
   return(K);
}


#' Function to compute a polynomial  kernel from a feature matrix 
#' @param W : a feature matrix, Rows are examples and columns are features      
#' @param degree : integer corresponding to a degree of the polynomial (def.: 2)
#' @param scale : scaling factor. Double larger than 0. If scale=-1 (def) scale is set to 1/ncol(W);
#' @param v : constant factor (def. 0)
#' @return A polynomial (REX:it was gaussian) kernel matrix 
poly.kernel <- function(W, degree=2, scale=-1, v=0) {
   n <- nrow(W);
   m <- ncol(W);
   nn <- rownames(W);
   if (scale == -1)
     scale <- 1/as.double(m);
   K <- matrix(numeric(n*n), nrow=n);
   K <- .C("poly_kernel", as.double(K), as.double(W), as.integer(degree), as.double(scale), as.double(v), as.integer(n), as.integer(m), PACKAGE = "RANKS")[[1]];
   K <- matrix(K, nrow=n);
   rownames(K) <- colnames(K) <- nn;
   return(K);
}





################################################################
## Score functions 
## N.B. Per questi metodi i primi 2 argomenti risultano invertiti rispetto alla precedente versione
################################################################

# eav.score  Empirical Average Score 

########################################################
#' Method to compute the Empirical Average score for a single vertex
#' \eqn{score(x) = - K(x,x) + 2/|x.pos| * sum_{x_i \in x.pos} K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' return The eav score of the element x
#' @export
#' @docType methods
#' @rdname single.eav.score-methods
setGeneric("single.eav.score", 
                 function(RW, x, x.pos, auto=FALSE) standardGeneric("single.eav.score"));

#' Method to compute the Empirical Average score for a single vertex
#' \eqn{score(x) = - K(x,x) + 2/|x.pos| * sum_{x_i \in x.pos} K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' return The eav score of the element x
setMethod("single.eav.score", signature(RW="matrix"),
  function(RW, x, x.pos, auto=FALSE) {
   if (nrow(RW) != ncol(RW))
     stop("second arg must be a square matrix");
   if (auto)
     score <- -RW[x,x] + (2/length(x.pos)) * sum(RW[x,x.pos])   
   else
     score <- (2/length(x.pos)) * sum(RW[x,x.pos]); 
   return(score);
})


#' Method to compute the Empirical Average score for a single vertex
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @return The eav score of the element x
setMethod("single.eav.score", signature(RW="graph"),
  function(RW, x, x.pos, auto=FALSE) {
      RW <- as(RW, "matrix");
      return(single.eav.score(RW, x, x.pos, auto));
})


########################################################
#' Method to compute the Empirical Average score for a set of vertices
#' \eqn{score(x) = - K(x,x) + 2/|x.pos| * sum_{x_i \in x.pos} K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return Vector of the eav scores of the elements x. The names of the vector correspond to the indices x
#' @export
#' @docType methods
#' @rdname eav.score-methods
setGeneric("eav.score", 
                 function(RW, x, x.pos, auto=FALSE, norm=TRUE) standardGeneric("eav.score"));

#' Method to compute the Empirical Average score for a set of vertices
#' \eqn{score(x) = - K(x,x) + 2/|x.pos| * sum_{x_i \in x.pos} K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return Vector of the eav scores of the elements x. The names of the vector correspond to the indices x
setMethod("eav.score", signature(RW="matrix"),
  function(RW, x, x.pos, auto=FALSE, norm=TRUE) {
	if (nrow(RW) != ncol(RW))
      stop("second arg must be a square matrix");
	if (auto)
      score <- -diag(RW)[x] + (2/length(x.pos)) * apply(as.matrix(RW[x,x.pos]),1,sum) 
	else
      score <- (2/length(x.pos)) * apply(as.matrix(RW[x,x.pos]),1,sum);
	names(score)<-x;  
	if (norm) {
	  ma <- max(score);
	  if (ma>0)
	    score <- score/ma;
	}
	return(score);
})


#' Method to compute the Empirical Average score for a set of vertices
#' @param x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed
#' @param  RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @return Vector of the eav scores of the elements x. The names of the vector correspond to the indices x
setMethod("eav.score", signature(RW="graph"),
  function(RW, x, x.pos, auto=FALSE, norm=TRUE) {
     RW <- as(RW, "matrix");
     return(eav.score(RW, x, x.pos, auto, norm));
})


#' Method to compute the Kernel NN score for a single vertex
#' \eqn{score(x) = - min_{x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#'              between pairs of nodes.
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @return The NN score of the element x
#' @export
#' @docType methods
#' @rdname single.NN.score-methods
setGeneric("single.NN.score", 
                 function(RW, x, x.pos, auto=FALSE) standardGeneric("single.NN.score"));
				 
#' Method to compute the Kernel NN score for a single vertex
#' \eqn{score(x) = - min_{x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @return The NN score of the element x
setMethod("single.NN.score", signature(RW="matrix"),
  function(RW, x, x.pos, auto=FALSE) {
	if (nrow(RW) != ncol(RW))
      stop("second arg must be a square matrix");
	n <- length(x.pos);
	if (auto) {
      scores <- rep(RW[x,x],n) + diag(RW)[x.pos] - 2 * RW[x,x.pos];
      names(scores) <- x;
      score <-  - min(scores);           
	} else {
      scores <- RW[x,x.pos];
      names(scores) <- x;
      score <- max(scores);
	}
	return(score);
})


#' Method to compute the Kernel NN score for a single vertex
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @return The NN score of the element x
setMethod("single.NN.score", signature(RW="graph"),
  function(RW, x, x.pos, auto=FALSE) {
     RW <- as(RW, "matrix");
     return(single.NN.score(RW, x, x.pos, auto));
})

#' Method to compute the Kernel NN score for a set of vertices
#' \eqn{score(x) = - min_{x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#' @param x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return The vector of the NN scores of the elements x. The names of the vector correspond to the indices x
#' @export
#' @docType methods
#' @rdname NN.score-methods
setGeneric("NN.score", 
                 function(RW, x, x.pos, auto=FALSE, norm=TRUE) standardGeneric("NN.score"));

#' Method to compute the Kernel NN score for a set of vertices
#' \eqn{score(x) = - min_{x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#'              between pairs of nodes.
#' @param x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return The vector of the NN scores of the elements x. The names of the vector correspond to the indices x
setMethod("NN.score", signature(RW="matrix"),
  function(RW, x, x.pos, auto=FALSE, norm=TRUE) {
   if (nrow(RW) != ncol(RW))
     stop("second arg must be a square matrix");
   n <- length(x.pos);
   if (auto) {
     scores <- matrix(rep(diag(RW)[x],n), ncol=n) +
     matrix(rep(diag(RW)[x.pos], length(x)), nrow=length(x), byrow=T) - 2 * RW[x,x.pos];
     score <- - apply(as.matrix(scores),1,min);
	 names(score)<-x;           
   } else {
     scores <- RW[x,x.pos];
     score <- apply(as.matrix(scores),1,max); 
	 names(score)<-x;
   }
   if (norm) {
	  ma <- max(score);
	  if (ma>0)
	    score <- score/ma;
   }
   return(score);
})

#' Method to compute the Kernel NN score for a set of vertices
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @return The vector of the NN scores of the elements x. The names of the vector correspond to the indices x
setMethod("NN.score", signature(RW="graph"),
  function(RW, x, x.pos, auto=FALSE, norm=TRUE) {
     RW <- as(RW, "matrix");
     return(NN.score(RW, x, x.pos, auto, norm));
})

#' Method to compute the Kernel KNN score for a single vertex
#' \eqn{score(x) = - sum_{k nearest x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#'              between pairs of nodes.
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @param k : integer. Number of the k nearest neighbours to be considered
#' @return The KNN score of the element x
#' @export
#' @docType methods
#' @rdname single.KNN.score-methods
setGeneric("single.KNN.score", 
                 function(RW, x, x.pos, k=3, auto=FALSE) standardGeneric("single.KNN.score"));

#' Method to compute the Kernel KNN score for a single vertex
#' \eqn{score(x) = - sum_{k nearest x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#'              between pairs of nodes.
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @param k : integer. Number of the k nearest neighbours to be considered
#' @return The KNN score of the element x
setMethod("single.KNN.score", signature(RW="matrix"),
  function(RW, x, x.pos, k=3, auto=FALSE) {
	if (nrow(RW) != ncol(RW))
      stop("second arg must be a square matrix");
	n <- length(x.pos);
	if (k > n) {
      k <- n;
      warn.message <- paste("k must be lower or equal to the number of positive examples: k set to", n);
      warning(warn.message);
	}     
	if (auto) {
      scores <- rep(RW[x,x],n) + diag(RW)[x.pos] - 2 * RW[x,x.pos];
      names(scores) <- x;
      scores <- sort(scores)[1:k];
      score <-  - sum(scores);           
	} else {
      scores <- RW[x,x.pos];
      names(scores) <- x;
      scores <- sort(scores, decreasing=TRUE)[1:k];
      score <- sum(scores);
	}
	return(score/k);
})

#' Method to compute the Kernel KNN score for a single vertex
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default)
#' @param k : integer. Number of the k nearest neighbours to be considered
#' @return The KNN score of the element x
setMethod("single.KNN.score", signature(RW="graph"),
  function(RW, x, x.pos, k=3, auto=FALSE) {
    RW <- as(RW, "matrix");
    return(single.KNN.score(RW, x, x.pos, k, auto));
})

#' Method to compute the Kernel KNN score for a set of vertices.
#' \eqn{score(x) = - (1/k) |sum_{k nearest x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#'              between pairs of nodes.
#' @param x : vector of integer. Indices correspond to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded
#' (default). The option TRUE is at the moment significantly slower.
#' @param k : integer. Number of the k nearest neighbours to be considered (def. 3)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return the KNN scores of the elements x
#' @export
#' @docType methods
#' @rdname KNN.score-methods
setGeneric("KNN.score", 
                 function(RW, x, x.pos, k=3, auto=FALSE, norm=TRUE) standardGeneric("KNN.score"));

#' Method to compute the Kernel KNN score for a set of vertices. 
#' \eqn{score(x) = - (1/k) |sum_{k nearest x_i \in x.pos} ( K(x,x) + K(x_i,x_i) -2 K(x,x_i)}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity
#'              between pairs of nodes.
#' @param x : vector of integer. Indices correspond to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded 
#' (default). The option TRUE is at the moment significantly slower.
#' @param k : integer. Number of the k nearest neighbours to be considered (def. 3)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return the KNN scores of the elements x
setMethod("KNN.score", signature(RW="matrix"),
  function(RW, x, x.pos, k=3, auto=FALSE, norm=TRUE) {
   #if (length(x.pos)==1)
    # stop("The number of positives must be larger than 1. This bug must be fixed");
   if (nrow(RW) != ncol(RW))
     stop("second arg must be a square matrix");
   # dyn.load("ranks.so");
   n <- length(x.pos);
   if (k==1)
     stop("k must be larger than 1. Please, use NN.score instead.");
   if (k > n) {
     k <- n;
     warn.message <- paste("k must be lower or equal to the number of positive examples: k set to", n);
     warning(warn.message);
   }     
   if (auto) {	 
	 scores <- matrix(rep(diag(RW)[x],n), ncol=n) +
     matrix(rep(diag(RW)[x.pos], length(x)), nrow=length(x), byrow=T) - 2 * RW[x,x.pos];
     score <- apply(scores,1,sort); # results are "per column"
     score <- score[1:k,]; 
     score <- - apply(score,2,sum);
     names(score) <- x;                   
   } else {
     scores <- as.matrix(RW[x,x.pos]);
	 y <- numeric(k);
	 score = matrix(numeric(length(x)*k), ncol=k);
	 for (j in 1:length(x))
        score[j,] <- .C("select_top", as.double(scores[j,]), as.double(y), as.integer(n), as.integer(k), PACKAGE = "RANKS")[[2]]; 
     score <- apply(score,1,sum);  # results are "per rows"
     names(score) <- x; 
   }
   score<-score/k;
   if (norm) {
	  ma <- max(score);
	  if (ma>0)
	    score <- score/ma;
   }
   return(score);
}) 


#' Method to compute the Kernel KNN score for a set of vertices. 
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : vector of integer. Indices correspond to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded 
#' (default). The option TRUE is at the moment significantly slower.
#' @param k : integer. Number of the k nearest neighbours to be considered (def. 3)
#' @return The KNN scores of the elements x
setMethod("KNN.score", signature(RW="graph"),
  function(RW, x, x.pos, k=3, auto=FALSE, norm=TRUE) {
    RW <- as(RW, "matrix");
    return(KNN.score(RW, x, x.pos, k, auto, norm));
}) 



##### WSLD score #################################
#' Method to compute the Weighted Sum with Linear Decay (WSLD) score for a single vertex
#' Let \eqn{K(x, x_jk)} be the kth rank order index w.r.t. \eqn{x_j \in V_C}, and \eqn{m=|V_C|}, then:
#' \eqn{score(x) = max_{x_i \in x_pos} K(x,x_i) + \sum_k=2^m [(1/(d * (k-1))) * K(x, x_jk)]}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @param d : integer. coefficient of linear decay (def. 2)
#' @return The WSLD score of the element x
#' @export
#' @docType methods
#' @rdname single.WSLD.score-methods
setGeneric("single.WSLD.score", 
                 function(RW, x, x.pos, d=2, auto=FALSE) standardGeneric("single.WSLD.score"));

#' Method to compute the Weighted Sum with Linear Decay (WSLD) score for a single vertex
#' Let \eqn{K(x, x_jk)} be the kth rank order index w.r.t. \eqn{x_j \in V_C}, and \eqn{m=|V_C|}, then:
#' \eqn{score(x) = max_{x_i \in x_pos} K(x,x_i) + \sum_k=2^m [(1/(d * (k-1))) * K(x, x_jk)]}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @param d : integer. coefficient of linear decay (def. 2)
#' @return The WSLD score of the element x
setMethod("single.WSLD.score", signature(RW="matrix"),
  function(RW, x, x.pos, d=2, auto=FALSE) {
   s <- double(1);
   if (nrow(RW) != ncol(RW))
     stop("first arg must be a square matrix");
   if (auto)
     score <- -RW[x,x] + .C("wsld2", as.double(s), as.double(RW[x,x.pos]), as.integer(length(x.pos)), as.integer(d), PACKAGE = "RANKS")[[1]]   
   else
     score <- .C("wsld2", as.double(s), as.double(RW[x,x.pos]), as.integer(length(x.pos)), as.integer(d), PACKAGE = "RANKS")[[1]]; 
   return(score);
})

#' Method to compute the Weighted Sum with Linear Decay (WSLD) score for a single vertex
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : integer. Index corresponding to the element of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the component \eqn{-K(x,x)} is computed, otherwise is discarded (default)
#' @param d : integer. coefficient of linear decay (def. 2)
#' @return The WSLD score of the element x
setMethod("single.WSLD.score", signature(RW="graph"),
  function(RW, x, x.pos, d=2, auto=FALSE) {
      RW <- as(RW, "matrix");
      return(single.WSLD.score(RW, x, x.pos, d, auto));
})



#' Method to compute the Weighted Sum with Linear Decay (WSLD) score for a set of vertices.
#' Let \eqn{K(x, x_jk)} be the kth rank order index w.r.t. \eqn{x_j \in V_C}, and \eqn{m=|V_C|}, then:
#' \eqn{score(x) = max_{x_i \in x_pos} K(x,x_i) + \sum_k=2^m [(1/(d * (k-1))) * K(x, x_jk)]}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : vector of integer. Indices correspond to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default). The option TRUE is at the moment significantly slower.
#' @param d : integer. coefficient of linear decay (def. 2)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return The WSLD scores of the elements x
#' @export
#' @docType methods
#' @rdname WSLD.score-methods
setGeneric("WSLD.score", 
                 function(RW, x, x.pos, d=2, auto=FALSE, norm=TRUE) standardGeneric("WSLD.score"));

#' Method to compute the Weighted Sum with Linear Decay (WSLD) score for a set of vertices. 
#' Let \eqn{K(x, x_jk)} be the kth rank order index w.r.t. \eqn{x_j \in V_C}, and \eqn{m=|V_C|}, then:
#' \eqn{score(x) = max_{x_i \in x_pos} K(x,x_i) + \sum_k=2^m [(1/(d * (k-1))) * K(x, x_jk)]}
#' @param RW : matrix. It must be a kernel matrix or a symmetric matrix expressing the similarity between pairs of nodes.
#' @param x : vector of integer. Indices correspond to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default). The option TRUE is at the moment significantly slower.
#' @param d : integer. coefficient of linear decay (def. 2)
#' @param norm : boolean. If TRUE (def.) the scores are normalized between 0 and 1.
#' @return The WSLD scores of the elements x
setMethod("WSLD.score", signature(RW="matrix"),
  function(RW, x, x.pos, d=2, auto=FALSE, norm=TRUE) {
  
    # Internal Function to compute the WSLD function (weighted sum with linear decay)
    # Input:
    # x : a numeric vector
    # d : the decay constant
    # N.B. : this version is not efficient. It must be implemented in C.
    wsld <- function (x,d) {
       ind.pos <- which(x>0);
       len <- length(ind.pos);
       if (len == 0)
    	 return (0);
       sorted <- sort(x[ind.pos], decreasing=TRUE);
       w <- d*(0:(len-1));
       w[1] <- 1;
       return(sum(sorted/w));
    }

   
    if (nrow(RW) != ncol(RW))
     stop("WSLD.score: first arg must be a square matrix");
	m <- length(x);
	score <- double(m);
   	if (auto)
      score <- -diag(RW)[x] + apply(as.matrix(RW[x,x.pos]),1,wsld,d=d) 
	else  {
	  M <- t(RW[x,x.pos]);
	  score <- .C("do_wsld_scores_from_matrix", as.double(score), as.double(M), as.integer(m), as.integer(length(x.pos)), as.integer(d), PACKAGE = "RANKS")[[1]];
	}
	names(score)<-x;  
	if (norm) {
	  ma <- max(score);
	  if (ma>0)
	    score <- score/ma;
    }
	return(score);
}) 


#' Method to compute the Weighted Sum with Linear Decay (WSLD) score for a set of vertices. 
#' @param RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
#' @param x : vector of integer. Indices correspond to the elements of the RW matrix for which the score must be computed
#' @param x.pos : vector of integer. Indices of the positive elements of the RW matrix
#' @param auto : boolean. If TRUE the components \eqn{K(x,x) + K(x_i,x_i)} are computed, otherwise are discarded (default). The option TRUE is at the moment significantly slower.
#' @param d : integer. coefficient of linear decay (def. 2)
#' @return The WSLD scores of the elements x
setMethod("WSLD.score", signature(RW="graph"),
  function(RW, x, x.pos, d=2, auto=FALSE, norm=TRUE) {
    RW <- as(RW, "matrix");
    return(WSLD.score(RW, x, x.pos, d, auto,norm));
}) 




#################################
## Functions for cross validation
#################################


######################################################################
## Method to perform a simple cross-validation with ranking

#' Function to perform a simple cross-validation with a kernel-based score method
#' @param RW : matrix. It could be a kerenl matrix or the adjacency matrix of a graph
#' @param ind.pos : indices of the positive examples. They are the indices the row of RW corresponding to positive examples.
#' @param m : number of folds (def: 5)
#' @param init.seed : initial seed for the random generator. If NULL (default) no initialization is performed
#' @param fun : function. It must be a kernel-based score method: KNN.score (default), NN.score or eav.score
#' @param ... : optional arguments for the function fun
#' @return a vector with the scores computed on each example
ker.score.cv <- function(RW, ind.pos, m=5, init.seed=NULL, fun=KNN.score, ...) {
   if (!is.null(init.seed))
     set.seed(init.seed);
   n <- nrow(RW);
   scores <- numeric(n);
   names(scores) <- rownames(RW);
   
   # Realization of the folds
   fold <- do.stratified.cv.data(1:n, ind.pos, k=m, seed=init.seed);
   
   # computing scores on the k folds
   for (i in 1:m) {
     x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
     core.pos <- integer(0);
     for (j in 1:m)
       if (j!=i)
	 core.pos <- c(core.pos, fold$fold.positives[[j]]);
     scores[x] <- fun(RW, x, core.pos, ...);
   }
   return(scores);
}


#####################################################################
#### Functions to perform multiple-cross validations and computing average scores
#####################################################################


#' Function to execute multiple cross-validation with a kernel-based score method for ranking
#' It computes the scores by averaging across multiple cross validations
#' @param RW : matrix. Random walk matrix or any valid symmetric matrix
#' @param ind.pos : indices of the positive examples. They are the indices the row of RW corresponding to positive examples.
#' @param m : number of folds for each cross-validation
#' @param p : number of repeated cross-validations
#' @param init.seed : initial seed for the random generator (def: 0)
#' @param fun : function. It must be a kernel-based score method: KNN.score (default), NN.score or eav.score
#' @param ... : optional arguments for the function fun
#' @return a list with two components:
#' - av.scores : a vector with the average scores across multiple cross-validations.
#'               Elements of the vector av.scores correspond to the rows of RW
#' - pos.scores : a vector with the scores of positive elements collected at each iteration
multiple.ker.score.cv <-  function(RW, ind.pos, m=5, p=100, init.seed=0, fun=KNN.score, ...) {
   n <- nrow(RW);
   current.scores <- av.scores <- numeric(n);  # vector of average scores
   pos.scores <- numeric();  # vector collecting the scores of positive examples
   
   for (v in 1:p) {
	 # Realization of the m folds
	 fold <- do.stratified.cv.data(1:n, ind.pos, k=m, seed=v+init.seed);
	 # computing scores on the m folds
	 for (i in 1:m) {
       x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
       core.pos <- integer(0);
       for (j in 1:m)
    	 if (j!=i)
	       core.pos <- c(core.pos, fold$fold.positives[[j]]);
       current.scores[x] <- fun(RW,x,core.pos, ...);	   
	 }
	 av.scores <- av.scores + current.scores;
	 pos.scores <- c(pos.scores,current.scores[ind.pos]);
   }
   return(list(av.scores=av.scores/p, pos.scores=pos.scores));
}





#' Function to find the optimal quantile alpha and corresponding threshold by  cross-validation with a kernel-based 
#' score method. The optimality is computed with respect to a specific metric( def: F-score).
#' @param K : matrix. Kernel matrix
#' @param ind.pos : indices of the positive examples. They are the indices the row of K corresponding to positive examples of the training set.
#' @param ind.non.pos : indices of the non positive examples. They are the indices the row of K corresponding to non positive examples  of the training set.
#' @param m : number of folds (default: 5)
#' @param alpha : vector of the quantiles to be tested
#' @param init.seed : initial seed for the random generator. If NULL (def) no initialization is performed
#' @param opt.fun : function. Function implementing the metric to choice the optimal threshold. The F-score (compute.F) is the default.
#'           Available functions:
#'                 - compute.F (default) : F-score
#'                 - compute.acc  : accuracy
#'           N.B.: any function having two arguments representing the vector of predicted and true labels
#'           can be in principle used.
#' @param fun : function. It must be a kernel-based score method:
#'                 - KNN.score (default)
#'                 - NN.score
#'                 - eav.score
#' @param ... : optional arguments for the function fun
#' @return a list with 3 elements:
#' alpha : quantile corresponding to the best F-score
#' thresh : threshold corresponding to the best F-score
#' pos.scores : scores of the positive elements computed through CV
find.optimal.thresh.cv <- function(K, ind.pos, ind.non.pos, m=5, 
   alpha=seq(from=0.05, to=0.6, by=0.05), init.seed=NULL, opt.fun=compute.F, fun=KNN.score, ...) {
   
   if (!is.null(init.seed))
     set.seed(init.seed);
   ind.all <- c(ind.pos, ind.non.pos);
   n=length(ind.all);
   K <- K[ind.all, ind.all];
   ind.pos <- 1:length(ind.pos);
   ind.non.pos <- (length(ind.pos)+1):n;
   av.scores <- numeric(n);  # vector of average scores
   pos.scores <- numeric();  # vector collecting the scores of positive examples cumulated across the folds
   fold <- do.stratified.cv.data(1:n, ind.pos, k=m, seed=init.seed);
   for (i in 1:m) {
     x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
     core.pos <- integer(0);
     for (j in 1:m)
       if (j!=i)
	     core.pos <- c(core.pos, fold$fold.positives[[j]]);
     av.scores[x] <- fun(K,x,core.pos, ...);	   
   }
   pos.scores <- av.scores[ind.pos];
   
   opt.F = -1;
   opt.alpha=0;
   labels <- c(rep(1, times=length(ind.pos)), rep(0, times=length(ind.non.pos)));
   for (i in 1:length(alpha)) {
     threshold <- selection.test(pos.scores, av.scores, ind.positives=ind.pos, alpha=alpha[i])$thresh;
     pred <- ifelse(av.scores>=threshold, 1, 0);
	 F <- opt.fun(pred,labels);    
	 if (F > opt.F) {
	   opt.F <- F;
	   opt.alpha <- alpha[i];
	   opt.thresh <- threshold;
	 }
   }
   return(list(alpha=opt.alpha, thresh=opt.thresh, pos.scores=pos.scores));
}



#' Function to execute multiple cross-validation with a kernel-based score method and to find the optimal
#' threshold for a given class by internal cross-validation. 
#' Scores are computed by averaging across multiple external cross-validations.
#' The optimal quantile and corresponding threshold  are selected by internal cross-validation using a 
#' specific metric (def: F-score).
#' @param K : matrix. Kernel matrix
#' @param ind.pos : indices of the positive examples. They are the indices the row of K corresponding to positive examples.
#' @param m : number of folds for each cross-validation (the number applies to both external and internal CV)
#' @param p : number of repeated external cross-validations
#' @param alpha : vector of the quantiles to be tested
#' @param init.seed : initial seed for the random generator (def: 0)
#' @param fun : function. It must be a kernel-based score method: KNN.score (default) or NN.score or eav.score
#' @param ... : optional arguments for the function fun
#' @return a list with three components:
#' - av.scores : a vector with the average scores across multiple cross-validations.
#'               Elements of the vector av.scores correspond to the rows of K
#' - opt.alpha : the optimal alpha for the class
#' - opt.thresh : the optimal threshold for the class
multiple.ker.score.thresh.cv <- function(K, ind.pos, m=5, p=100, alpha=seq(from=0.05, to=0.6, by=0.05), init.seed=0, fun=KNN.score, ...) {
   set.seed(init.seed);
   n <- nrow(K);
   current.scores <- av.scores <- numeric(n);  # vector of average scores
   opt.alpha = 0;
   opt.thresh = 0;
   for (v in 1:p) {
	 # Realization of the m folds
	 fold <- do.stratified.cv.data(1:n, ind.pos, k=m, seed=v+init.seed);
	 # computing scores on the m folds
	 for (i in 1:m) {
       x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
       core.pos <- core.non.pos <- integer(0);
       for (j in 1:m)
    	 if (j!=i) {
	       core.pos <- c(core.pos, fold$fold.positives[[j]]);
	       core.non.pos <- c(core.non.pos, fold$fold.non.positives[[j]]);
	   }
       current.scores[x] <- fun(K,x,core.pos, ...);
	   res <- find.optimal.thresh.cv(K, core.pos, core.non.pos, 
	                                m=m, alpha=alpha, init.seed=init.seed, fun=fun, ...);
	   opt.alpha <- opt.alpha + res$alpha;	
	   opt.thresh <- opt.thresh + res$thresh;   	   
	 }
	 av.scores <- av.scores + current.scores;
   }
   return(list(av.scores=av.scores/p, opt.alpha=opt.alpha/(m*p), opt.thresh=opt.thresh/(m*p)));
}


#' Function to classify labels according to an hold-out procedure with a kernel-based score method.
#' The optimal threshold for a given class is obtained by (possibly multiple) internal cross-validation. 
#' Scores of the held-out nodes are computed and thresholds computed on the training set 
#' by cross-validation are then used to classify the held-out nodes
#' The optimal quantile and corresponding threshold  are selected by internal cross-validation using the F-score as metrics
#' @param K : matrix. Kernel matrix
#' @param ind.pos : indices of the positive examples of the training set. They are the indices the row of K corresponding to
#'          positive examples of the training set
#' @param ind.test : indices of the examples of the test set. They are the indices the row of K corresponding to
#'           examples of the test set
#' @param m : number of folds for each cross-validation 
#' @param p : number of repeated  cross-validations on the training set
#' @param alpha : vector of the quantiles to be tested
#' @param init.seed : initial seed for the random generator (def: 0)
#' @param opt.fun : function. Function implementing the metric to choice the optimal threshold.
#'           The F-score (compute.F) is the default.
#'           Available functions:
#'                 - compute.F (default) : F-score
#'                 - compute.acc  : accuracy
#'           N.B.: any function having two arguments representing the vector of predicted and true labels
#'           can be in principle used.
#' @param fun : function. It must be a kernel-based score method:
#'                 - KNN.score (default)
#'                 - NN.score 
#'                 - eav.score
#' @param ... : optional arguments for the function fun
#' @return a list with four components:
#' - labels : vector of the predicted labels for the test set(1 represent positive, 0 negative)
#' - scores : a vector with the  scores computed on the test set.
#'               Elements of the vector av.scores correspond to ind.test rows of K
#' - opt.alpha : the optimal alpha for the class
#' - opt.thresh : the optimal threshold for the class
ker.score.classifier.holdout <- function(K, ind.pos, ind.test, m=5, p=10, alpha=seq(from=0.05, to=0.6, by=0.05), init.seed=0, opt.fun=compute.F, fun=KNN.score, ...) {
   set.seed(init.seed);
   n <- nrow(K);
   ind.train <- (1:n)[-ind.test];
   dd <- intersect(ind.pos, ind.test);
   if (length(dd)>0)
     stop("ker.score.classifier.holdout: conflicting indices between training and test set");
   opt.alpha = 0;
   opt.thresh = 0;
   for (v in 1:p) {
	 # Realization of the m folds
	 fold <- do.stratified.cv.data(ind.train, ind.pos, k=m, seed=v+init.seed);
	 # computing thresholds on the m folds
	 for (i in 1:m) {
       x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
       core.pos <- core.non.pos <- integer(0);
       for (j in 1:m)
    	 if (j!=i) {
	       core.pos <- c(core.pos, fold$fold.positives[[j]]);
	       core.non.pos <- c(core.non.pos, fold$fold.non.positives[[j]]);
	   }
	   res <- find.optimal.thresh.cv(K, core.pos, core.non.pos, 
	                                m=m, alpha=alpha, init.seed=init.seed, opt.fun=opt.fun, fun=fun, ...);
	   opt.alpha <- opt.alpha + res$alpha;	
	   opt.thresh <- opt.thresh + res$thresh;   	   
	 }
   }
   opt.alpha <- opt.alpha/(m*p); 
   opt.thresh <- opt.thresh/(m*p);
   scores <- fun(K,ind.test, ind.pos, ...);
   labels <- labelsfromscores(scores, opt.thresh);
   names(labels) <- names(scores) <- names(ind.test);
   return(list(labels=labels, scores=scores, opt.alpha=opt.alpha, opt.thresh=opt.thresh));
}




#' Function to classify labels accroding to external cross-validation procedure with a kernel-based score method and to find the optimal
#' threshold for a given class by internal cross-validation. 
#' Scores are computed by averaging across (possibly) multiple external cross-validations.
#' The optimal quantile and corresponding threshold  are selected by internal cross-validation using the F-score as metrics
#' @param K : matrix. A Kernel matrix
#' @param ind.pos : indices of the positive examples. They are the indices the row of K corresponding to positive examples.
#' @param m : number of folds for each cross-validation (the number applies to both external and internal CV)
#' @param p : number of repeated external cross-validations
#' @param alpha : vector of the quantiles to be tested
#' @param init.seed : initial seed for the random generator (def: 0)
#' @param opt.fun : function. Function implementing the metric to choice the optimal threshold.
#'           The F-score (compute.F) is the default.
#'           Available functions:
#'                 - compute.F (default) : F-score
#'                 - compute.acc  : accuracy
#'           N.B.: any function having two arguments representing the vector of predicted and true labels
#'           can be in principle used.
#' @param fun : function. It must be a kernel-based score method:
#'                 - KNN.score (default)
#'                 - NN.score 
#'                 - eav.score
#' @param ... : optional arguments for the function fun
#' @return a list with four components:
#' - labels : vector of the predicted labels (1 represent positive, 0 negative)
#' - av.scores : a vector with the average scores across multiple cross-validations.
#'               Elements of the vector av.scores correspond to the rows of K
#' - opt.alpha : the optimal alpha for the class
#' - opt.thresh : the optimal threshold for the class
ker.score.classifier.cv <- function(K, ind.pos, m=5, p=100, alpha=seq(from=0.05, to=0.6, by=0.05), init.seed=0, opt.fun=compute.F, fun=KNN.score, ...) {
   set.seed(init.seed);
   n <- nrow(K);
   current.scores <- av.scores <- numeric(n);  # vector of average scores
   opt.alpha = 0;
   opt.thresh = 0;
   for (v in 1:p) {
	 # Realization of the m folds
	 fold <- do.stratified.cv.data(n, ind.pos, k=m, seed=v+init.seed);
	 # computing scores on the m folds
	 for (i in 1:m) {
       x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
       core.pos <- core.non.pos <- integer(0);
       for (j in 1:m)
    	 if (j!=i) {
	       core.pos <- c(core.pos, fold$fold.positives[[j]]);
	       core.non.pos <- c(core.non.pos, fold$fold.non.positives[[j]]);
	   }
       current.scores[x] <- fun(K, x, core.pos, ...);
	   res <- find.optimal.thresh.cv(K, core.pos, core.non.pos, 
	                                m=m, alpha=alpha, init.seed=init.seed, opt.fun=opt.fun, fun=fun, ...);
	   opt.alpha <- opt.alpha + res$alpha;	
	   opt.thresh <- opt.thresh + res$thresh;   	   
	 }
	 av.scores <- av.scores + current.scores;
   }
   av.scores <- av.scores/p;
   opt.alpha <- opt.alpha/(m*p); 
   opt.thresh <- opt.thresh/(m*p);
   labels <- labelsfromscores(av.scores, opt.thresh)
   return(list(labels=labels, av.scores=av.scores, opt.alpha=opt.alpha, opt.thresh=opt.thresh));
}





###########################################################################
## Non parameteric selection test based on multiple cross-validations
###########################################################################

#' Non parametric test to select the most significant unlabeled examples.
#' @param pos.scores : vector with scores of positive examples. It is returned from multiple.ker.score.cv.
#' @param av.scores : a vector with the average scores computed by multiple.ker.score.cv. It may be a named vector. If not, the names attributes corresponding to the indices of the vector are added.
#' @param ind.positives : indices of the positive examples. They are the indices of av.scores corresponding to positive examples.
#' @param alpha : significance level (def. 0.05)
#' @param thresh.pos : only values larger than thresh.pos are retained in pos.scores (def.: 0)
#' @return a list with 5 components:
#' - selected : a named vector with the components of av.scores selected by the test
#' - selected.labeled : a named vector with the labeled components of av.scores selected by the test
#' - selected.unlabeled : a named vector with the unlabeled components of av.scores selected by the test
#' - thresh : the score threshold selected by the test
#' - alpha : significance level (the same value of the input)
selection.test <- function(pos.scores, av.scores, ind.positives, alpha=0.05, thresh.pos=0) {
  
  n <- length(av.scores);
  pos.scores <- pos.scores[pos.scores > thresh.pos];
  if (is.null(names(av.scores)))
    names(av.scores) <- 1:n; 
  av.scores.labeled <- av.scores[ind.positives];
  av.scores.unlabeled <- av.scores[-ind.positives];
  thresh <- quantile(pos.scores,alpha);
  selected <- av.scores[av.scores>=thresh];
  selected.labeled <- av.scores.labeled[av.scores.labeled>=thresh];
  selected.unlabeled <- av.scores.unlabeled[av.scores.unlabeled>=thresh];
  
  return(list(selected=selected, selected.labeled=selected.labeled,
              selected.unlabeled=selected.unlabeled, thresh=thresh, alpha=alpha));
}


############# Utility functions ############################# 


######################################################################
#' Function to generate data for the stratified cross-validation.
#' @param examples : indices of the examples (a vector of integer)
#' @param positives : vector of integer. Indices of the positive examples. The indices refer to the indices of examples
#' @param k : number of folds (def = 5)
#' @param seed : seed of the random generator (def=NULL). If is set to NULL no initiazitation is performed
#' @return a list with 2 two components
#'   - fold.non.positives : a list with k components. Each component is a vector with the indices of the non positive elements of the fold
#'   - fold.positives : a list with k components. Each component is a vector with the indices of the positive elements of the fold
#' N.B.: in both elements indices refer to row numbers of the data matrix	 
do.stratified.cv.data <- function(examples, positives, k=5, seed=NULL) {
  if (!is.null(seed))
     set.seed(seed);
  fold.non.positives <- fold.positives <- list();
  for (i in 1:k) {
    fold.non.positives[[i]] <- integer(0);
    fold.positives[[i]] <- integer(0);
  }
  # examples <- 1:n;
  non.positives <- setdiff(examples,positives);
  # non.positives <- examples[-positives];
  non.positives <- sample(non.positives);
  positives <- sample(positives);
  n.positives <- length(positives);
  resto.positives <- n.positives%%k;
  n.pos.per.fold <- (n.positives - resto.positives)/k;
  n.non.positives <- length(non.positives);
  resto.non.positives <- n.non.positives%%k;
  n.non.pos.per.fold <- (n.non.positives - resto.non.positives)/k;
  j=1; 
  if (n.non.pos.per.fold > 0)
    for (i in 1:k) {
      fold.non.positives[[i]] <- non.positives[j:(j+n.non.pos.per.fold-1)];
      j <- j + n.non.pos.per.fold;
    }
  j.pos=1;  
  if (n.pos.per.fold > 0)
    for (i in 1:k) {
      fold.positives[[i]] <- positives[j.pos:(j.pos+n.pos.per.fold-1)];
      j.pos <- j.pos + n.pos.per.fold;
    }
  
  if (resto.non.positives > 0)
    for (i in k:(k-resto.non.positives+1)) {
      fold.non.positives[[i]] <- c(fold.non.positives[[i]], non.positives[j]);
      j <- j + 1;
    }
  
  if (resto.positives > 0) 
    for (i in 1:resto.positives) {
      fold.positives[[i]] <- c(fold.positives[[i]], positives[j.pos]);
      j.pos <- j.pos + 1;
    }
  
  return(list(fold.non.positives=fold.non.positives, fold.positives=fold.positives));
}




#' Function to compute the F-measure for a single class
#' @param pred : factor of the predicted labels
#' @param labels : factor of the true labels
#' Note that 0 level stands for negative and 1 for positive.
#' In general the first level is negative and the second positive
#' @return The computed F-score
compute.F <- function(pred,labels) {      
     if (length(pred)!=length(labels))
         stop("compute.F: lengths of true and predicted labels do not match.");
     neg.ex <- which(labels <= 0);	
	 np <- length(neg.ex);
	 pos.ex <- which(labels > 0);
	 npos <- length(pos.ex);	 
	 TP <- sum(pred[pos.ex] > 0);
	 FN <- sum(pred[pos.ex] <= 0);	
	 TN <- sum(pred[neg.ex] <= 0);
	 FP <- sum(pred[neg.ex] > 0);	           
     if ((TP+FP) == 0)
       precision <- 0
     else 
       precision <- TP/(TP+FP);
     if ((TP+FN) == 0)
       recall <- 0
     else
       recall <- TP/(TP+FN);
     if ((precision+recall) == 0)
     	F <- 0
     else
     	F = 2 *(precision*recall) / (precision+recall);      
     return (F);
}

#' Function to compute the accuracy for a single class
#' @param pred : vector of the predicted labels
#' @param labels : vector of the true labels
#' Note that 0  stands for negative and 1 for positive.
#' In general the first level is negative and the second positive
#' @return The computed accuracy
compute.acc <- function(pred,labels) {      
     if (length(pred)!=length(labels))
         stop("compute.acc: lengths of true and predicted labels do not match.");
     n <- length(labels);
	 acc <- sum(pred==labels)/n;
	 return(acc);
}



#' Function to compute the labels of multiple classes from the corresponding scores
#' @param S : matrix. Matrix of scores: rows represent examples, columns classes
#' @param thresh : numeric vector. Vector of the thesholds for each class.
#' @return a binary matrix with the labels of the predictions: rows represent examples, columns classes. Element L[i,j] is the label of example i w.r.t. class j.  L[i,j]=1 if i belongs to j, 0 otherwise.
Multiple.labels.from.scores <- function(S, thresh.vect) {
   n <- length(thresh.vect);
   L <- matrix(integer(nrow(S)*ncol(S)), nrow=nrow(S));
   for (i in 1:n) 
     L[,i] <- ifelse(S[,i]>=thresh.vect[i], 1, 0);
   rownames(L) <- rownames(S);
   colnames(L) <- colnames(S);
   return(L);
}

#' Function to compute the labels of a single class from the corresponding scores
#' @param scores : numeric. Vector of scores: each element correspond to the score of an example
#' @param thresh : real value. Threshold for the class.
#' @return Labels of the predictions: a vector with 0 or 1 value. The label res[i]=1 if scores[i]>thresh, otherwise res[i]=0
labelsfromscores <- function(scores, thresh) {
   res <- ifelse(scores>=thresh, 1, 0);
   names(res) <- names(scores);
   return(res);
}

#' Function that computes the norm 1 of a numeric vector
#' @param x : numeric vector
#' @return a single real value (the norm1 of the input vector)
norm1 <- function(x) {
 return(sum(abs(x)));
}

#' Function to normalize a kernel according to the unit sphere
#' @param K : a kernel matrix     
#' @return The normalized kernel
Unit.sphere.norm <- function(K) {
   d <- diag(K);
   return(K / sqrt(d %*% t(d)));
}



######################################################################################################################
# HIGH LEVEL FUNCTIONS
#######################################################################################################################

#' do.RANKS
#' High level function to perform experiments with RANKS
#' It perform a k fold CV repeated rep times on a given data set
#' @param score : function. It must be a kernel-based score method:
#'                 - eav.score (default)
#'                 - NN.score 
#'                 - KNN.score
#'                 - WSLD.score
#' @param kernel : kernel method (def. rw.kernel)
#' @param a : kernel parameter (def. 2)
#' @param p : number of steps of the RW kernel (def. 1)
#' @param sparsify : boolean. If TRUE (def) the input matrix is sparsified using Sparsify.matrix from NetpreProc
#' @param kk : number of folds of the cross validation (def: 5)
#' @param rep : number of repetitions of the cross validation (def: 1)
#' @param seed : intialization seed for the random generator to create folds (def:0)
#' @param data.dir : relative path to directory where the adjiacency matrix is stored (def: ../data)
#' @param labels.dir : relative path to directory where the label matrix is stored (def: ../data)
#' @param output.dir : relative path to directory where the results are stored  (def: ../Results)
#' @param data : name of the data set to loaded (without rda extension). It must be  an .rda file containing the adjiacency matrix of the graph. 
#'        It assumes that it is in the data.dir directory
#' @param labels : name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples.
#'          It assumes that it is in the labels.dir directory. Note that data and labels must have the same number of rows and in the same order
#' @param ... : optional arguments to be passed to the fucntion multiple.ker.score.cv that performs the CV
#' @return 3 rda files stored in the Results directory:
#' - Scores results: A matrix with examples on rows and classes on columns representing the computed scores for each example and for each considered class
#' - AUC results files computed through AUC.single.over.classes
#' - Precision at given recall results computed through precision.at.multiple.recall.level.over.classes
do.RANKS  <- function(score=eav.score, kernel=rw.kernel, a=2,  p=1,  sparsify=TRUE, kk=5, rep=1, seed=0, 
                       data.dir="../data/", labels.dir="../data/", output.dir="../Results/", data, labels, ...)  {
  
  recall.levels <- c(0.01, 0.05, seq(from=0.1, to=1, by=0.1));
  
  # loading the adjacency matrix
  dataset.name <- paste0(data.dir,data,".rda");
  data.name=load(dataset.name);
  K=eval(parse(text=data.name));  # K represents the adjacency matrix
  K[K<0]<-0;

  # loading labels matrix
  dataset.name <- paste0(labels.dir,labels,".rda");
  label.name=load(dataset.name);
  T=eval(parse(text=label.name));  # T represents the label matrix
  classnames<-colnames(T);
  if ("GO:0008150" %in% classnames) {
  T <- T[,-which(classnames=="GO:0008150")]
  } else if ("GO:0003674" %in% classnames) {
  T <- T[,-which(classnames=="GO:0003674")]
  } else if ("GO:0005575" %in% classnames) {
  T <- T[,-which(classnames=="GO:0005575")]
  } else if ("HP:0000001" %in% classnames) {
  T <- T[,-which(classnames=="HP:0000001")]
  };

  nclasses <- ncol(T);
  nelem <- nrow(T);
  k=d=0; # to avoid warnings from check --as-cran
  
  if (nelem != nrow(K))
    stop("do.RANKS: adjacency and label matrices do not agree");
  
  if (sparsify)
    K <- Sparsify.matrix(K, 1);
  K <- kernel(K, a);
  if (p>1)
    K <- p.step.rw.kernel(K, p=p);

  # construction of the matrix of scores. 
  S  <- matrix(numeric(nclasses * nelem), nrow=nelem);
  elemnames<-rownames(T);
  classnames<-colnames(T);
  rownames(S) <- elemnames;
  colnames(S) <- classnames;

  # computing scores for each class
  
  
  
  for (i in 1:nclasses)  {
      ind.pos <- which(T[,i]==1);
      res <- multiple.ker.score.cv(K, ind.pos, m=kk, p=rep, init.seed=seed, fun=score, ...)$av.scores; 
      S[,i] <- res;
      cat("Class ", i, " : ", classnames[i], "\n");
  }
  
  # saving scores
  score.name <- as.character(substitute(score));
  arg.dots = list(...);
  if ((score.name == "KNN.score") && hasArg(k))
    score.name <- paste0(score.name, arg.dots$k)
  else if (score.name == "KNN.score")
    score.name <- paste0(score.name, 3)
  else if ((score.name == "WSLD.score") && hasArg(d))
    score.name <- paste0(score.name, arg.dots$d)
  else if (score.name == "WSLD.score")
	score.name <- paste0(score.name, 2);
	
  score.file = paste0(output.dir, "Scores.",score.name,".","p",p,".","a",a,".","f",sparsify,".",data.name, ".", label.name, ".rda");
  save(S, file=score.file)
  
  # computing and saving AUC  
  
  AUC <- AUC.single.over.classes(T, S);
  
  AUC.file = paste0(output.dir, "AUC.",score.name,".","p",p,".","a",a,".","f",sparsify,".",data.name, ".", label.name, ".rda");
  
  save(AUC, file=AUC.file);
  
  # computing and saving PXR 
  PXR <- precision.at.multiple.recall.level.over.classes (T, S, rec.levels=recall.levels);
  
  PXR.file = paste0(output.dir, "PXR.",score.name,".","p",p,".","a",a,".","f",sparsify,".",data.name, ".", label.name, ".rda");
  
  save(PXR, file=PXR.file);
  
}


#' do.loo.RANKS
#' High level function to perform leave one out (loo) experiments with RANKS
#' It perform a loo on a given data set
#' @param score : function. It must be a kernel-based score method:
#'                 - KNN.score 
#'                 - NN.score 
#'                 - eav.score  (default)
#'                 - WSLD
#' @param compute.kernel : logical. If TRUE (def.) a kernel matrix is computed from data according to the choice of the function kernel, otherwise the data matrix is used as it is.
#' @param kernel : kernel method (def. rw.kernel)
#' @param a : kernel parameter (def. 2)
#' @param k : number of neighbours for KNN.score. It is meaningful only for  kNN  (def.19)
#' @param d : integer. coefficient of linear decay for the WSLD score. It is meaningful only for  the WSLD score  (def.2) 
#' @param p : number of steps of the RW kernel (def. 1)
#' @param sparsify : boolean. If TRUE the input matrix is sparsified using Sparsify.matrix from NetpreProc (def: FALSE)
#' @param norm : logical. If TRUE for each class the score is normalized in [0,1], otherwise the raw scores are maintained (default).
#' @param data : name of the network data set to be loaded (without rda extension). It must be  an .rda file containing the adjiacency matrix of the graph. 
#'        By default it assumes that it is in the "data" directory
#' @param labels : name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples.
#'          By default it assumes that it is in the "data" directory
#' @param output.name : name of the output file (without rda extension). Other informations including the learning parameters are added
#' @param net.dir : relative path to directory where the adjiacency matrix is stored (def: data)
#' @param labels.dir : relative path to directory where the label matrix is stored (def: data)
#' @param output.dir : relative path to directory where the results are stored  (def: Results). Note that data and labels must have the same number of rows and in the same order. Moreover if any label column corresponds to any GO root term, this is eliminated to avoid prediction of GO root nodes.
#' @return 3 rda files stored in the Results directory:
#' - Scores results: A matrix with examples on rows and classes on columns representing the computed scores for each example and for each considered class
#' - AUC results files computed through AUC.single.over.classes
#' - Precision at given recall results computed through precision.at.multiple.recall.level.over.classes
do.loo.RANKS  <- function(score=eav.score, compute.kernel=TRUE, kernel=rw.kernel,  a=2, k=19, d=2, 
                 p=1, sparsify=FALSE, norm=FALSE, data, labels, output.name, 
				 net.dir="data/", labels.dir="data/", output.dir="Results/")  {

  recall.levels <- c(0.01, 0.05, seq(from=0.1, to=1, by=0.1));
  
  # loading the adjacency matrix
  dataset.name <- paste0(net.dir, data,".rda");
  data.name=load(dataset.name);
  M=eval(parse(text=data.name));  # M represents the adjacency matrix
  M[M<0]<-0;

  # loading labels matrix
  dataset.name <- paste0(labels.dir,labels,".rda");
  label.name=load(dataset.name);
  ann=eval(parse(text=label.name));  # ann represents the label matrix
  # if there is a GO root node or a HPO root node this is eliminated
  classnames<-colnames(ann);
  if ("GO:0008150" %in% classnames) {
  ann <- ann[,-which(classnames=="GO:0008150")]
  } else if ("GO:0003674" %in% classnames) {
  ann <- ann[,-which(classnames=="GO:0003674")]
  } else if ("GO:0005575" %in% classnames) {
  ann <- ann[,-which(classnames=="GO:0005575")]
  } else if ("HP:0000001" %in% classnames) {
  ann <- ann[,-which(classnames=="HP:0000001")]
  };


  nclasses <- ncol(ann);
  nelem <- nrow(ann);
  
  if (any(rownames(M) != rownames(ann)))
    stop("do.loo.RANKS: adjacency and label matrices do not agree");
  
  if (sparsify)
     M <- Sparsify.matrix(M, k=1);
  if (compute.kernel) {
    M <- kernel(M, a);
    if (p>1)
      M <- p.step.rw.kernel(M, p=p);
  }
  diag(M) <- 0;  # necessary for loo, otherwise AUC=1

  # construction of the matrix of scores. 
  S  <- matrix(numeric(nclasses * nelem), nrow=nelem);
  elemnames<-rownames(ann);
  classnames<-colnames(ann);
  rownames(S) <- elemnames;
  colnames(S) <- classnames;

  # computing scores for each class
  score.name <- as.character(substitute(score));
  
  if (score.name == "KNN.score")
    for (i in 1:nclasses)  {
      ind.pos <- which(ann[,i]==1);
	  res <- KNN.score(M, 1:nelem, ind.pos, k=k, auto=FALSE, norm=norm);
      S[,i] <- res;
      cat("Class ", i, " : ", classnames[i], "\n");
    }
  else if (score.name == "WSLD.score")
    for (i in 1:nclasses)  {
      ind.pos <- which(ann[,i]==1);
	  res <- WSLD.score(M, 1:nelem, ind.pos, d=d, auto=FALSE, norm=norm);
      S[,i] <- res;
      cat("Class ", i, " : ", classnames[i], "\n");
    }
  else 
    for (i in 1:nclasses)  {
      ind.pos <- which(ann[,i]==1);
	  res <- score(M, 1:nelem, ind.pos, auto=FALSE, norm=norm);
      S[,i] <- res;
      cat("Class ", i, " : ", classnames[i], "\n");
    }
  
  rm(M);gc();
  # saving scores
  if (score.name == "KNN.score")
    score.name <- paste0(score.name, k);
  
  
  score.file = paste0(output.dir,"Scores.",score.name,".","p",p,".","a",a,".",output.name, ".rda");
  save(S, file=score.file)
  
  # computing and saving AUC  
  
  AUC <- AUC.single.over.classes(ann, S);
  
  AUC.file = paste0(output.dir,"AUC.loo.",score.name,".","p",p,".","a",a,".",output.name, ".rda");
  
  save(AUC, file=AUC.file);
  
  # computing and saving PXR 
  PXR <- precision.at.multiple.recall.level.over.classes (ann, S, rec.levels=recall.levels);
  
  PXR.file = paste0(output.dir,"PXR.loo.",score.name,".","p",p,".","a",a,".",output.name, ".rda");
  
  save(PXR, file=PXR.file);
  rm(ann,S);gc();  
}





######################################################################
# Implementation of Random walk, Random walk with restart, labelprop and GBA methods
######################################################################


#' Function that performs a random Walk on a given graph
#' @param W : adjacency matrix of the graph
#' @param ind.positives : indices of the "core" positive examples of the graph.
#'                They represent to the indices of W corresponding to the positive examples
#' @param tmax : maximum number of iterations (def: 1000)
#' @param eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-10)
#' @param norm : if TRUE (def) the adjacency matrix W of the graph is normalized to M = D^-1 * W, otherwise
#'        it is assumed that the matrix W is just normalized
#' @return a list with three elements:
#' - p : the probability at the steady state or after tmax iterations
#' - ind.positives: indices of the "core" positive examples of the graph (it is equal to the same
#'                  input parameter
#' - n.iter : number of performed iterations
RW <- function(W, ind.positives, tmax=1000, eps=1e-10, norm=TRUE) {
   if (norm) 
     M <-Prob.norm(W) # M = D^-1 * W
   else
     M <- W;
   n <- nrow(M);
   p0 <- p <- numeric(n);
   names(p) <- names(p0) <- rownames(W);
   rm(W);
   n.positives <- length(ind.positives);
   if (n.positives == 0)
 	  stop("RW: number of core positives is equal to 0!");
   p0[ind.positives] <- 1/n.positives;
   
   #cat("Inizio trasposta \n"); print(date());
   #M <- t(M);
   #cat("Fine trasposta \n"); print(date());
   p <- p0;
   for (t in 1:tmax) {
 	 p0 <- p;
 	 p <- as.vector(p0 %*% M);   # no explicit transpose is computed: this is equivalent to t(M) %*% p0
 	 if (norm1(p-p0) < eps) break(); 
   }  
   return(list(p=p, ind.positives=ind.positives, n.iter=t));  
}





#' Function that performs a random Walk with restart (RWR) on a given graph
#' @param W : adjacency matrix of the graph
#' @param ind.positives : indices of the "core" positive examples of the graph. They represent to the indices of W corresponding to the positive examples
#' @param gamma : restart parameter (def: 0.6)
#' @param tmax : maximum number of iterations (def: 1000)
#' @param eps : maximum allowed difference between the computed probabilities at the steady state
#' @param norm : if TRUE (def) the adjacency matrix W of the graph is normalized to M = D^-1 * W, otherwise it is assumed that the matrix W is just normalized
#' @return a list with three elements:
#' - p : the probability at the steady state
#' - ind.positives : indices of the "core" positive examples of the graph (it is equal to the same
#'                  input parameter
#' - n.iter : number of performed iterations
RWR <- function(W, ind.positives, gamma=0.6, tmax=1000, eps=1e-10, norm=TRUE) {  
   if (norm) 
     M <-Prob.norm(W) # M = D^-1 * W
   else
     M <- W;
   n <- nrow(M);
   p0 <- p <- numeric(n);
   names(p) <- names(p0) <- rownames(W);
   rm(W);
   n.positives <- length(ind.positives);
   if (n.positives == 0)
 	  stop("RWR: number of core positives is equal to 0!");
   p0[ind.positives] <- 1/n.positives;
   
   # M <- t(M);
   p <- p0;
   for (t in 1:tmax) {
 	 pold <- p;
 	 p <- ((1-gamma) * as.vector(pold %*% M)) + gamma * p0;
 	 if (norm1(p-pold) < eps) break(); 
   }  
   return(list(p=p, ind.positives=ind.positives, n.iter=t));  
}


#' Function that implements the Label propagation algorithm of Zhu and Ghahramani
#' @param W : adjacency matrix of the graph
#' @param ind.positives : indices of the "core" positive examples of the graph.
#'                They represent to the indices of W corresponding to the positive examples
#' @param tmax : maximum number of iterations (def: 1000)
#' @param eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-5)
#' @param norm : if TRUE (def) the adjacency matrix W of the graph is normalized to M = D^-1 * W, otherwise
#'        it is assumed that the matrix W is just normalized
#' @return a list with three elements:
#' - p : the probability at the steady state
#' - ind.positives: indices of the "core" positive examples of the graph (it is equal to the same
#'                  input parameter)
#' - n.iter : number of performed iterations
label.prop <- function(W, ind.positives, tmax=1000, eps=1e-5, norm=TRUE) { 
   if (norm) 
     M <-Prob.norm(W) # M = D^-1 * W
   else
     M <- W;
   n <- nrow(M);
   p <- numeric(n);
   names(p)  <- rownames(W);
   n.positives <- length(ind.positives);
   if (n.positives == 0)
 	  stop("label.prop: number of core positives is equal to 0!");
   p[ind.positives] <- 1;
   
   M <- t(M);
   for (t in 1:tmax) {
 	 pold <- p;
 	 p <- M %*% pold;
	 p[ind.positives] <- 1;
 	 if (norm1(p-pold) < eps) break();  
   }  
   return(list(p=p, ind.positives=ind.positives, n.iter=t));  
}

#' Function that implements a simple Guilt By Association (GBA) method for  label ranking based on 
#' the sum of edge weights connecting a node to its positive neighbours
#' @param W : adjacency matrix of the graph
#' @param ind.positives : indices of the "core" positive examples of the graph. They represent to the indices of W corresponding to the positive examples.
#' @return a list with one element:
#' - p : score associated to each node
GBAsum <- function(W, ind.positives) {
   n <- nrow(W);
   p <- numeric(n);
   names(p)  <- rownames(W);
   n.positives <- length(ind.positives);
   if (n.positives == 0)
 	  stop("GBAsum: number of core positives is equal to 0!");
   diag(W)<-0;
   p <- apply(W[,ind.positives],1,sum);
   res <- list(p=p);
   return(res);
}

#' Function that implements a simple Guilt By Association (GBA) method for  label ranking based on 
#' the maximum between the edge weights connecting a node to its positive neighbours
#' @param W : adjacency matrix of the graph
#' @param ind.positives : indices of the "core" positive examples of the graph. They represent to the indices of W corresponding to the positive examples.
#' @return a list with one element :
#' - p (score associated to each node)
GBAmax <- function(W, ind.positives) {
   n <- nrow(W);
   p <- numeric(n);
   names(p)  <- rownames(W);
   n.positives <- length(ind.positives);
   if (n.positives == 0)
 	  stop("GBAmax: number of core positives is equal to 0!");
   diag(W)<-0;
   p <- apply(W[,ind.positives],1,max);
   res <- list(p=p);
   return(res);
}


######################################################################
## Method to perform a simple cross-validation 

#' Function to execute cross-validation with random walk based, labelprop and GBA methods
#' @param W : adjacency matrix of the graph
#'     note that if the optional argument norm=TRUE (def.), the W matrix is normalized, otherwise it
#'     is assumed that W is just normalized
#' @param ind.pos : indices of the "core" positive examples of the graph. They represent to the indices of W corresponding to the positive examples
#' @param k : number of folds (def: 5)
#' @param init.seed : initial seed for the random generator. If 0 (default) no initialization is performed
#' @param fun : function. It must be a randow walk method:
#'                 - RW (default)
#'                 - RWR
#'                 - label.prop
#'                 - GBAsum
#'                 - GBAmax
#' @param ... : optional arguments for the function fun:
#' - gamma : restart parameter (def: 0.6) (meaningful only for RWR)
#' - tmax : maximum number of iterations (def: 1000)
#' - eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-10)
#' @return a vector with the the probabilities for each example at the steady state

RW.cv <- function(W, ind.pos, k=5, init.seed=0, fun=RW, ...) {
   if (init.seed!=0)
     set.seed(init.seed);
   n <- nrow(W);
   p <- numeric(n);
   names(p) <- rownames(W);
   
   # Realization of the folds. ERA QUI L'ERRORE: prima al posto di 1:n c'era solo n !!!!
   fold <- do.stratified.cv.data(1:n, ind.pos, k=k, seed=init.seed);
   
   # computing scores on the k folds
   for (i in 1:k) {
     x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
     core.pos <- integer(0);
     for (j in 1:k)
       if (j!=i)
	     core.pos <- c(core.pos, fold$fold.positives[[j]]);
     p[x] <- fun(W,core.pos, ...)$p[x];
   }
   return(p);
}


######################################################################
## Method to perform multiple cross-validation

#' Function to execute multiple cross-validation with random walk based, labelprop and GBA methods.
#' It computes the scores by averaging across multiple cross validations
#' @param W : adjacency matrix of the graph
#'     note that if the optional argument norm=TRUE (def.), the W matrix is normalized, otherwise it
#'     is assumed that W is just normalized
#' @param ind.pos : indices of the "core" positive examples of the graph.
#'                They represent to the indices of W corresponding to the positive examples
#' @param k : number of folds (def: 5)
#' @param p : number of repeated cross-validations
#' @param init.seed : initial seed for the random generator. If 0 (default) no initialization is performed
#' @param fun : function. It must be a randow walk method:
#'                 - RW (default)
#'                 - RWR
#'                 - label.prop
#'                 - GBAsum
#'                 - GBAmax
#' @param ... : optional arguments for the function fun:
#' - gamma : restart parameter (def: 0.6) (meaningful only for RWR)
#' - tmax : maximum number of iterations (def: 1000)
#' - eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-10)
#' @return a vector with the the probabilities for each example at the steady state averaged across multiple cross-validations
multiple.RW.cv <- function(W, ind.pos, k=5, p=100, init.seed=0, fun=RW, ...) {
   n <- nrow(W);
   current.p <- average.p <- numeric(n);  # vector of average scores
   
   for (v in 1:p) {
	 # Realization of the k folds
	 fold <- do.stratified.cv.data(1:n, ind.pos, k=k, seed=v+init.seed);
	 # computing scores on the k folds
	 for (i in 1:k) {
       x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]]);
       core.pos <- integer(0);
       for (j in 1:k)
    	 if (j!=i)
	   core.pos <- c(core.pos, fold$fold.positives[[j]]);
       current.p[x] <- fun(W,core.pos, ...)$p[x]; 
	 }
	 average.p <- average.p + current.p;
   }
   return(average.p/p);
}

#######################################################################
# High level functions to perform experiments with  RW, RWR  and GBA
#######################################################################

#' do.RWR
#' High level function to perform experiments with RWR
#' It perform a k fold CV repeated 1 time on a given data set
#' @param gamma : restart parameter (def: 0.6) 
#' @param tmax : maximum number of iterations (def: 1000)
#' @param eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-10)
#' @param k : number of folds for the cross validation (def. 5)
#' @param filter : if TRUE (def) the adjacency matrix is sparsified otherwise not
#' @param seed : seed of the random generator for the generation of the folds (def: 1):
#' @param data : name of the data set to loaded (without rda extension). It must be  an .rda file containing the adjiacency matrix of the graph. 
#'        It assumes that it is in the "data" directory
#' @param labels : name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples.
#'          It assumes that it is in the "data" directory
#' @return 3 rda files stored in the Results directory (names of the files are automatically computed from the input data):
#' - Scores results: A matrix with examples on rows and classes on columns representing the computed scores for each example and for each considered class
#' - AUC results files computed through AUC.single.over.classes
#' - Precision at given recall results computed through precision.at.multiple.recall.level.over.classes
do.RWR  <- function(gamma=0.6, tmax=1000, eps=1e-10, k=5, filter=TRUE, seed=1, data, labels)  {

  recall.levels <- c(0.01, 0.05, seq(from=0.1, to=1, by=0.1));
  
  # loading the adjacency matrix
  dataset.name <- paste0("data/",data,".rda");
  data.name=load(dataset.name);
  K=eval(parse(text=data.name));  # K represents the adjacency matrix
  K[K<0]<-0;

  # loading labels matrix
  dataset.name <- paste0("data/",labels,".rda");
  label.name=load(dataset.name);
  T=eval(parse(text=label.name));  # T represents the label matrix

  nclasses <- ncol(T);
  ndrugs <- nrow(T);
  
  if (ndrugs != nrow(K))
    stop("do.RWR: adjacency and label matrices do not agree");
  
  if (filter)
    K <- Sparsify.matrix(K, k=1); 
  K <-Prob.norm(K); 

  # construction of the matrix of scores. 
  S  <- matrix(numeric(nclasses * ndrugs), nrow=ndrugs);
  drugnames<-rownames(T);
  classnames<-colnames(T);
  rownames(S) <- drugnames;
  colnames(S) <- classnames;

  # computing scores for each class
  for (i in 1:nclasses)  {
    ind.pos <- which(T[,i]==1);
    # 1 CV
    res <- RW.cv(K, ind.pos, k=k, init.seed=seed, fun=RWR, gamma=gamma, tmax=tmax, eps=eps, norm=FALSE); 
    S[,i] <- res;
    cat("Class ", i, " : ", classnames[i], "\n");
  }

  # saving scores
  
  score.file = paste0("Results/Scores.","RWR.","g",gamma,".","f",filter,".",data.name, ".", label.name, ".rda");
  save(S, file=score.file);
  
  # computing and saving AUC  
  
  AUC <- AUC.single.over.classes(T, S);
  
  AUC.file = paste0("Results/AUC.","RWR.","g",gamma,".","f",filter,".",data.name, ".", label.name, ".rda");  
  save(AUC, file=AUC.file);

  # computing and saving PXR 
  PXR <- precision.at.multiple.recall.level.over.classes (T, S, rec.levels=recall.levels);

  PXR.file = paste0("Results/PXR.","RWR.","g",gamma,".","f",filter,".",data.name, ".", label.name, ".rda");  
  
  save(PXR, file=PXR.file);
  
}


#' do.RW
#' High level function to perform experiments with RW
#' It perform a k fold CV repeated 1 time on a given data set
#' @param tmax : maximum number of iterations (def: 1000)
#' @param eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-10)
#' @param k : number of folds for the cross validation (def. 5)
#' @param filter : if TRUE (def) the adjacnecy matrix is sparsified otherwise not
#' @param seed : seed of the random generator for the generation of the folds (def: 1):
#' @param data : name of the data set to loaded (without rda extension). It must be  an .rda file containing the adjiacency matrix of the graph. 
#'        It assumes that it is in the "data" directory
#' @param labels : name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples.
#'          It assumes that it is in the "data" directory
#' @return 3 rda files stored in the Results directory (names of the files are automatically computed from the input data):
#' - Scores results: A matrix with examples on rows and classes on columns representing the computed scores for each example and for each considered class
#' - AUC results files computed through AUC.single.over.classes
#' - Precision at given recall results computed through precision.at.multiple.recall.level.over.classes
do.RW  <- function(tmax=1000, eps=1e-10, k=5, filter=TRUE, seed=1, data, labels)  {

  recall.levels <- c(0.01, 0.05, seq(from=0.1, to=1, by=0.1));
  
  # loading the adjacency matrix
  dataset.name <- paste0("data/",data,".rda");
  data.name=load(dataset.name);
  K=eval(parse(text=data.name));  # K represents the adjacency matrix
  K[K<0]<-0;

  # loading labels matrix
  dataset.name <- paste0("data/",labels,".rda");
  label.name=load(dataset.name);
  T=eval(parse(text=label.name));  # T represents the label matrix
  nclasses <- ncol(T);
  ndrugs <- nrow(T);
  
  if (ndrugs != nrow(K))
    stop("do.RW: adjacency and label matrices do not agree");
  
  if (filter)
    K <- Sparsify.matrix(K, k=1); 
  K <-Prob.norm(K); 

  # construction of the matrix of scores. 
  S  <- matrix(numeric(nclasses * ndrugs), nrow=ndrugs);
  drugnames<-rownames(T);
  classnames<-colnames(T);
  rownames(S) <- drugnames;
  colnames(S) <- classnames;

  # computing scores for each class
  for (i in 1:nclasses)  {
    ind.pos <- which(T[,i]==1);
    # 1 CV
    res <- RW.cv(K, ind.pos, k=k, init.seed=seed, fun=RW, tmax=tmax, eps=eps, norm=FALSE); 
    S[,i] <- res;
    cat("Class ", i, " : ", classnames[i], "\n");
  }

  # saving scores
  
  score.file = paste0("Results/Scores.","RW.", tmax, "step.","f",filter,".",data.name, ".", label.name, ".rda");
  save(S, file=score.file);
  
  # computing and saving AUC  
  
  AUC <- AUC.single.over.classes(T, S);
  
  AUC.file = paste0("Results/AUC.","RW.", tmax, "step.","f",filter,".",data.name, ".", label.name, ".rda");
  save(AUC, file=AUC.file);

  # computing and saving PXR 
  PXR <- precision.at.multiple.recall.level.over.classes (T, S, rec.levels=recall.levels);

  PXR.file = paste0("Results/PXR.","RW.", tmax, "step.","f",filter,".",data.name, ".", label.name, ".rda");
  
  save(PXR, file=PXR.file);
  
}

#' do.GBA
#' High level function to perform experiments with GBA
#' It perform a k fold CV repeated 1 time on a given data set
#' @param fun : function performing GBA. it can be one of the following
#'       - GBAsum: it sums the edge weights connecting a node to its positive neighbours
#'       - GBAmax: it computes the maximum between the edge weights connecting a node to its positive neighbours
#' @param k : number of folds for the cross validation (def. 5)
#' @param filter : if TRUE (def) the adjacnecy matrix is sparsified otherwise not
#' @param seed : seed of the random generator for the generation of the folds (def: 1):
#' @param data : name of the data set to loaded (without rda extension). It must be  an .rda file containing the adjiacency matrix of the graph. 
#'        It assumes that it is in the "data" directory
#' @param labels : name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples.
#'          It assumes that it is in the "data" directory
#' @return 3 rda files stored in the Results directory (names of the files are automatically computed from the input data):
#' - Scores results: A matrix with examples on rows and classes on columns representing the computed scores for each example and for each considered class
#' - AUC results files computed through AUC.single.over.classes
#' - Precision at given recall results computed through precision.at.multiple.recall.level.over.classes
do.GBA  <- function(fun=GBAsum, k=5, filter=TRUE, seed=1, data, labels)  {

  recall.levels <- c(0.01, 0.05, seq(from=0.1, to=1, by=0.1));
  
  # loading the adjacency matrix
  dataset.name <- paste0("data/",data,".rda");
  data.name=load(dataset.name);
  K=eval(parse(text=data.name));  # K represents the adjacency matrix
  K[K<0]<-0;
  gc();
  # loading labels matrix
  dataset.name <- paste0("data/",labels,".rda");
  label.name=load(dataset.name);
  T=eval(parse(text=label.name));  # T represents the label matrix
  gc();
  nclasses <- ncol(T);
  ndrugs <- nrow(T);
  
  if (ndrugs != nrow(K))
    stop("do.GBA: adjacency and label matrices do not agree");
  
  if (filter)
    K <- Sparsify.matrix(K, k=1); 
  gc();

  # construction of the matrix of scores. 
  S  <- matrix(numeric(nclasses * ndrugs), nrow=ndrugs);
  drugnames<-rownames(T);
  classnames<-colnames(T);
  rownames(S) <- drugnames;
  colnames(S) <- classnames;
  gc();
  
  # computing scores for each class
  for (i in 1:nclasses)  {
    ind.pos <- which(T[,i]==1);
    # 1 CV
    res <- RW.cv(K, ind.pos, k=k, init.seed=seed, fun=fun); 
    S[,i] <- res;
    cat("Class ", i, " : ", classnames[i], "\n");
	gc();
  }

  # saving scores
  
  fun.name <- as.character(substitute(fun));
  score.file = paste0("Results/Scores.", fun.name, ".f", filter, ".",data.name, ".", label.name, ".rda");
  save(S, file=score.file);
  
  # computing and saving AUC  
  
  AUC <- AUC.single.over.classes(T, S);
  
  AUC.file = paste0("Results/AUC.", fun.name, ".f", filter, ".",data.name, ".", label.name, ".rda");
  save(AUC, file=AUC.file)

  # computing and saving PXR 
  PXR <- precision.at.multiple.recall.level.over.classes (T, S, rec.levels=recall.levels);

  PXR.file = paste0("Results/PXR.", fun.name, ".f", filter, ".",data.name, ".", label.name, ".rda");
  
  save(PXR, file=PXR.file);
  
}






  
.onLoad <- function(libname=.libPaths(), pkgname="RANKS")
       library.dynam("RANKS", pkgname, libname);

.onAttach <- function(libname=.libPaths(), pkgname="RANKS")
               packageStartupMessage("RANKS: ranking and classification with kernelized score functions.\n")
