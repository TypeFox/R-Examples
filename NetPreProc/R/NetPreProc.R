# Methods to pre-process and normalize networked data
# May 2011
# October 2011
# December 2011

# dyn.load("NetPreProc.so");
# indx<<-0;

########################################################
setGeneric("Magnify.binary.features.norm", 
                 function(M) standardGeneric("Magnify.binary.features.norm"));

# Normalization of binary matrices according to Mostafavi et al. 2008
# Having a binary matrix M, for each feature, if b is the proportion of 1, then ones are replaced with -log(b) and zeros with log(1-b).
# Rows are examples and columns features to be normalized.
# Input:
# M : input binary matrix. Rows are examples, columns features
# Output :
# the pre-processed matrix
setMethod("Magnify.binary.features.norm", signature(M="matrix"),
  function(M) {
	n.feat <- ncol(M);
	n.ex <- as.double(nrow(M));
	M.mo <- matrix(numeric(nrow(M)*ncol(M)),nrow=nrow(M));
	for (i in 1:n.feat) {
	  b <- sum(M[,i])/n.ex;
	  M.mo[M[,i]==1, i] <- -log(b);
	  M.mo[M[,i]==0, i] <- log(1-b);	
	}
	# setting to  positive all weights
	M.mo <- M.mo - min(M.mo);
	rownames(M.mo) <- rownames(M);
	return(M.mo);
  }
)




########################################################
setGeneric("Chua.norm", 
                 function(W) standardGeneric("Chua.norm"));

# Normalization according to Chua et al., 2007
# The normalized weigth W_{ij} between nodes i and j is computed by taking into account their neighborhods N_i and N_j :
# W_{ij} = \frac{2|N_i \cap N_j|}{|N_i \setminus N_j| + 2|N_i \cap N_j| + 1}\times \frac{2|N_i \cap N_j|}{|N_j \setminus N_i| + 2|N_i \cap N_j| + 1}
# where $N_k$ is the set of the neighbors of gene $k$ ($k$ is included). 
# It is meaningful in particular with interaction data
# Input:
# W : input weight square symmetrix matrix of double
# Output :
# the pre-processed matrix
setMethod("Chua.norm", signature(W="matrix"),
function(W) {
  n = nrow(W);
  nn <- rownames(W);
  W <- .C("chua_like_norm", as.double(W), as.integer(n), PACKAGE="NetPreProc")[[1]];
  W <- matrix(W, ncol=n);
  rownames(W) <- colnames(W) <- nn;
  return(W);
})

# Normalization according to Chua et al., 2007
# The normalized weigth W_{ij} between nodes i and j is computed by taking into account their neighborhods N_i and N_j :
# W_{ij} = \frac{2|N_i \cap N_j|}{|N_i \setminus N_j| + 2|N_i \cap N_j| + 1}\times \frac{2|N_i \cap N_j|}{|N_j \setminus N_i| + 2|N_i \cap N_j| + 1}
# where $N_k$ is the set of the neighbors of gene $k$ ($k$ is included). 
# It is meaningful in particular with interaction data
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph.
# Output :
# the pre-processed matrix
setMethod("Chua.norm", signature(W="graph"),
function(W) {
  W <- as(W, "matrix");
  return(Chua.norm(W));
})


########################################################
setGeneric("Laplacian.norm", 
                 function(W) standardGeneric("Laplacian.norm"));

# Graph Laplacian normalization.
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by dividing each entry W_ij by the square root of the product of the sum of elements of row i and the sum of the elemnts in column j
# In other words if D is a diagonal matrix such that D_ii = \sum_j W_ij then the normalize matrix is: W_norm = D^(-1/2) * W * D^(-1/2)
# Input:
# W : symmetric matrix 
# Output:
# a normalized symmetric matrix
setMethod("Laplacian.norm", signature(W="matrix"),
function(W) {
# computing D^(-1/2) * W * D^(-1/2)
  n <- nrow(W);
  if (n != ncol(W))
    stop("first arg must be a square matrix");
  name.examples <- rownames(W);   
  diag.D <- apply(W,1,sum);
  diag.D[diag.D==0] <- Inf;
  inv.sqrt.diag.D <- 1/sqrt(diag.D);
  W <- .C("norm_lapl_graph", as.double(W), as.double(inv.sqrt.diag.D), as.integer(n), PACKAGE="NetPreProc")[[1]];
  W <- matrix(W, nrow=n);
  rownames(W) <- colnames(W) <- name.examples;
  return(W);
})


# Graph Laplacian normalization.
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by dividing each entry W_ij by the square root of the product of the sum of elements of row i and the sum of the elemnts in column j
# In other words if D is a diagonal matrix such that D_ii = \sum_j W_ij then the normalize matrix is: W_norm = D^(-1/2) * W * D^(-1/2)
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# Output:
# a normalized symmetric matrix
setMethod("Laplacian.norm", signature(W="graph"),
function(W) {
# computing D^(-1/2) * W * D^(-1/2)
  W <- as(W, "matrix");
  return(Laplacian.norm(W));
})



########################################################
setGeneric("Prob.norm", 
                 function(W) standardGeneric("Prob.norm"));

# Probabilistic normalization of a graph.
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by dividing each entry W_ij by the the sum of elements of row i 
# In other words if D is a diagonal matrix such that D_ii = \sum_j W_ij then the normalize matrix is: W_norm = D^(-1) * W 
# Input:
# W : symmetric matrix 
# Output:
# a normalized  matrix
# N.B.: La matrice risultante non e' simmetrica!
setMethod("Prob.norm", signature(W="matrix"),
function(W) {
   n <- nrow(W);
   if (n != ncol(W))
     stop("first arg must be a square matrix");
   names.var <- rownames(W);
   diag.D <- apply(W,1,sum);
   diag.D[diag.D==0] <- Inf;
   inv.diag.D <- 1/diag.D;
   W <- .C("norm_2", as.double(W), as.double(inv.diag.D), as.integer(n), PACKAGE="NetPreProc")[[1]];
   W <- matrix(W, nrow=n);
   rownames(W) <- colnames(W) <- names.var;
   return(W);
})

# Probabilistic normalization of a graph.
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by dividing each entry W_ij by the the sum of elements of row i 
# In other words if D is a diagonal matrix such that D_ii = \sum_j W_ij then the normalize matrix is: W_norm = D^(-1) * W 
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# Output:
# a normalized  matrix
# N.B.: La matrice risultante non e' simmetrica!
setMethod("Prob.norm", signature(W="graph"),
function(W) {
  W <- as(W, "matrix");
  return(Prob.norm(W));
})


########################################################
setGeneric("Max.Min.norm", 
                 function(W) standardGeneric("Max.Min.norm"));

# Graph normalization with respect to the minimum and maximum value of its weights
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by normalized each entry in the range [0..1]
# W.norm = (W - min(W))/(max(W)-min(W)) 
# Input:
# W : symmetric matrix 
# Output:
# a normalized symmetric matrix
setMethod("Max.Min.norm", signature(W="matrix"),
function(W) {
   min.W <- min(W);
   max.W <- max(W);
   return ((W-min.W)/(max.W-min.W));
})

# Graph normalization with respect to the minimum and maximum value of its weights
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by normalized each entry in the range [0..1]
# W.norm = (W - min(W))/(max(W)-min(W)) 
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# Output:
# a normalized symmetric matrix
setMethod("Max.Min.norm", signature(W="graph"),
function(W) {
   W <- as(W, "matrix");
   return(Max.Min.norm(W));
})


########################################################
# Function to obtain the Pearson correlation matrix between rows of a given matrix.
# You can also "sparsify" the matrix, by putting to 0 all the weights, by setting a threshold
# such that at least one edge is maintained for each node
# N.B.: The diagonal values are set to 0
# Input:
# M : input matrix
# cut : if TRUE (def.) at least one edge is maintained for each node, all the other edges are set to 0. If false
#       no threshold is computed
# remove.negatives : if TRUE (def) negative values are replaced with 0 in the correlation matrix
# min.thresh: minimum allowed threshold (def. 0). If a threshold lower than min.thresh is selected, than
#             it is substituted by min.thresh. 
#             Warning: setting min.thresh to high values may lead to highly disconneted network
# Output:
# a square symmetric matrix of the Pearson correlation coefficients 
# computed between the rows of M
Do.sim.matrix.Pearson <- function(M, cut=TRUE, remove.negatives=TRUE, min.thresh=0) {
   M.pc <- cor(t(M));
   M.pc[is.na(M.pc)|is.nan(M.pc)]<-0;
   diag(M.pc) <- 0;
   if (cut) {
     x <- apply(M.pc,1,max);  
     x <- x[x>0];  # remove 0 from maxima
     thresh <- min(x);
     if (thresh < min.thresh)
       thresh <- min.thresh;
     M.pc[M.pc<thresh] <- 0;   
   }
   if (remove.negatives && any(M.pc<0))
     M.pc[M.pc < 0] <- 0;
   rownames(M.pc) <- colnames(M.pc) <- rownames(M);
   return(M.pc);
}


########################################################
setGeneric("Sparsify.matrix", 
                 function(W, k=1) standardGeneric("Sparsify.matrix"));

# Method to sparsify a network matrix.
# By this function you can select a minimum of k edges for each node
# Input:
# W : a network matrix
# k : the number of guaranteed edges for each node (def.=1)
# Output:
# a square symmetric sparsified matrix 
setMethod("Sparsify.matrix", signature(W="matrix"),
function(W, k=1) {
   diag(W) <- 0;
   if (k==1) {
     x <- apply(W,1,max);
	 thresh <- min(x);
   } else {
     x <- apply(W,1,sort,TRUE);
     thresh <- min(x[k,]);
   }
   W[W<thresh] <- 0;  
   return(W);
})


# Method to sparsify a network matrix.
# By this function you can select a minimum of k edges for each node
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# k : the number of guaranteed edges for each node (def.=1)
# Output:
# a square symmetric sparsified matrix 
setMethod("Sparsify.matrix", signature(W="graph"),
function(W, k=1) {
   W <- as(W, "matrix");
   return(Sparsify.matrix(W,k));
})


########################################################
setGeneric("Sparsify.matrix.fixed.neighbours", 
                 function(W, k=10) standardGeneric("Sparsify.matrix.fixed.neighbours"));

# Method to sparsify a network matrix by fixing the number of edges for each node
# It selects the first k neighbours for each node (by row) according to the weight of the edge
# By this function you select exactly k edges for each node (if there are at least k edges in the adjacency matrix)
# Input:
# W : a network matrix
# k : the number of  edges for each node (def.=10)
# Output:
# a sparsified matrix (Warning: the matrix is not symmetric)
setMethod("Sparsify.matrix.fixed.neighbours", signature(W="matrix"),
function(W, k=10) {
   # diag(W) <- 0;
   n <- nrow(W);
   W <- apply(W,1, function(x) {
     x[rank(x, ties.method="max")<=(n-k)] <- 0;
     return(x);   
   });
   return(t(W));
})

# Method to sparsify a network matrix by fixing the number of edges for each node
# It selects the first k neighbours for each node (by row) according to the weight of the edge
# By this function you can select exactly k edges for each node (if there are at least k edges in the adjacency matrix)
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# k : the number of  edges for each node (def.=10)
# Output:
# a sparsified matrix (Warning: the matrix is not symmetric)
setMethod("Sparsify.matrix.fixed.neighbours", signature(W="graph"),
function(W, k=10) {
   W <- as(W, "matrix");
   return(Sparsify.matrix.fixed.neighbours(W,k));
})


########################################################
setGeneric("Binary.matrix.by.thresh", 
                 function(W, thresh=0.5) standardGeneric("Binary.matrix.by.thresh"));

# Method to transform a network matrix into a binary matrix.
# Values above the given threshold are set to 1, otherwise to 0
# Input:
# W : a network matrix
# thresh : the threshold (def.=0.5)
# Output:
# a binary matrix 
setMethod("Binary.matrix.by.thresh", signature(W="matrix"),
function(W, thresh=0.5) {
   diag(W) <- 0;
   W[W<=thresh] <- 0; 
   W[W>thresh] <- 1; 
   return(W);
})

# Method to transform a network matrix into a binary matrix.
# Values above the given threshold are set to 1, otherwise to 0
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# thresh : the threshold (def.=0.5)
# Output:
# a binary matrix 
setMethod("Binary.matrix.by.thresh", signature(W="graph"),
function(W, thresh=0.5) {
   W <- as(W, "matrix");
   return(Binary.matrix.by.thresh(W,thresh));
})


########################################################
setGeneric("check.network", 
                 function(W, name="Network matrix") standardGeneric("check.network"));
				 
# Method to check the characteristics of a network matrix
# Input:
# W : a network matrix
# name : a character vector that will be printed as heading
setMethod("check.network", signature(W="matrix"),
function(W, name="Network matrix") {
   cat(name, ": *** \n");
   if (any(is.na(W)))
     cat("WARNING: the matrix contains NA  ****", "\n");
   if (any(is.nan(W)))
     cat("WARNING:", name, "contains NaN  ****", "\n");
   if (any(is.infinite(W)))
     cat("WARNING: the matrix contains Inf values  ****", "\n");
   rr <- range(W);
   if ((rr[2] - rr[1]) <= 0)
     cat("WARNING: the matrix has values in a wrong range: ", rr, "***** \n");
   if (!isSymmetric(W))
     cat("WARNING: the matrix is not symmetric  ****** \n");
   dimen <- dim(W);
   edges <- sum(W>0);
   perc.edges <- round((edges/(dimen[1]*dimen[2])),5);
   cat("the network matrix has dimension: ", dimen, "\n");
   cat("the network has edge values in the range: ", rr, "\n");
   cat("the network has number of edges: ", edges, "\n");
   cat("the network has density: ", perc.edges, "\n");
})


# Method to check the characteristics of a network matrix
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# name : a character vector that will be printed as heading
setMethod("check.network", signature(W="graph"),
function(W, name="Network matrix") {
   W <- as(W, "matrix");
   return(check.network(W,name));
})

.onLoad <- function(libname=.libPaths(), pkgname="NetPreProc")
	   library.dynam("NetPreProc", pkgname, libname);
      
.onAttach <- function(libname=.libPaths(), pkgname="NetPreProc")
              packageStartupMessage("NetPreProc: Network Pre-Processing package for graph normalization. \n");
			  
