#' Convert weights, as COO matrix, to CSC matrix and weights vector
#'
#' \code{convert} Takes a weights triple (COO matrix) and converts it to a sparse edge-incidence matrix (CSC matrix) and weights vector.
#'
#' @param W COO matrix of weights: (i,j,w[ij])
#' @import Matrix
#' @export
#' @examples
#' W <- matrix(0,3,3)
#' W[1,] <- c(1,2,1)
#' W[2,] <- c(1,3,2)
#' W[3,] <- c(2,3,3)
#'
#' sol <- convert(W)
convert <- function(W) {
  m <- nrow(W)
  i <- j <- x <- integer(2*m)
  i[1:m] <- i[(m+1):(2*m)] <- 1:m
  j[1:m] <- W[,1]
  j[(m+1):(2*m)] <- W[,2]
  x[1:m] <- rep(1,m)
  x[(m+1):(2*m)] <- rep(-1,m)
  Phi <- sparseMatrix(i=i, j=j, x=x)
  w <- W[,3]
  return(list(Phi=Phi, w=w))
}

## Clusterpath preprocessing
tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}

vec2tri <- function(k,n) {
  i <- ceiling(0.5*(2*n-1 - sqrt((2*n-1)^2 - 8*k)))
  j <- k - n*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
}

#' Construct indices matrices
#'
#' \code{compactify_edges} constructs M1, M2, and ix index matrices.
#' @param w weights vector
#' @param n number of points to cluster
compactify_edges <- function(w,n) {
  sizes1 <- double(n)
  sizes2 <- double(n)

  P <- vec2tri(w@i+1,n)
  nEdge <- nrow(P)

  M1 <- matrix(0,nEdge,n)
  M2 <- matrix(0,nEdge,n)

  for (i in 1:n) {
    group1 <- which(P[,1] == i)
    sizes1[i] <- length(group1)
    if (sizes1[i] > 0) {
      M1[1:sizes1[i],i] <- group1
    }
    group2 <- which(P[,2] == i)
    sizes2[i] <- length(group2)
    if (sizes2[i] > 0) {
      M2[1:sizes2[i],i] <- group2
    }
  }

  M1 <- M1[1:max(sizes1),,drop=FALSE]
  M2 <- M2[1:max(sizes2),,drop=FALSE]

  return(list(ix=P,M1=M1,M2=M2,s1=sizes1,s2=sizes2))
}

#' Gaussian Kernel + k-Nearest Neighbor Weights
#'
#' \code{gkn_weights} combines Gaussian kernel weights with k-nearest neighbor weights
#'
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param phi The nonnegative parameter that controls the scale of kernel weights
#' @param k_row The number of row nearest neighbors
#' @param k_col The number of column nearest neighbors
#' @export
gkn_weights <- function(X,phi=0.5,k_row=5,k_col=5) {
  p <- nrow(X); n <- ncol(X)
  ## Construct Gaussian kernel weights
  w_row <- kernel_weights(t(X),phi/n)
  w_col <- kernel_weights(X,phi/p)
  ## Thin weights to k-nearest neighbors
  w_row <- knn_weights(w_row,k_row,p)
  w_col <- knn_weights(w_col,k_col,n)
  ## Normalize weights to sum to 1
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  ## Rescale weights to make column and row penalties commensurate
  w_row <- w_row/sqrt(n)
  w_col <- w_col/sqrt(p)

  ## Construct edge-incidence matrices
  E_row <- create_edge_incidence(w_row,p)
  E_col <- create_edge_incidence(w_col,n)

  ## Get connectivity information
  nRowComp <- length(find_clusters(weights_graph(w = w_row,p))$size)
  nColComp <- length(find_clusters(weights_graph(w = w_col,n))$size)

  return(list(w_row=w_row@x,w_col=w_col@x,E_row=E_row,E_col=E_col,
              nRowComp=nRowComp,nColComp=nColComp))
}

#' "Thin" a weight vector to be positive only for its k-nearest neighbors
#'
#' \code{knn_weights} takes a weight vector \code{w} and sets the ith
#' component \code{w[i]} to zero if either of the two corresponding nodes
#' is not among the other's \code{k} nearest neighbors.
#'
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param k The number of nearest neighbors
#' @param n The number of data points.
#' @return A vector \cite{w} of weights for convex clustering.
knn_weights <- function(w,k,n) {
  i <- 1
  neighbors <- tri2vec(i,(i+1):n,n)
  keep <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(n-1)) {
    group_A <- tri2vec(i,(i+1):n,n)
    group_B <- tri2vec(1:(i-1),i,n)
    neighbors <- c(group_A,group_B)
    knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep <- union(knn,keep)
  }
  i <- n
  neighbors <- tri2vec(1:(i-1),i,n)
  knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep <- union(knn,keep)
  if (length(keep) > 0)
    w[-keep] <- 0
  return(Matrix(data=w,ncol=1,sparse=TRUE))
}

#' Compute Gaussian Kernel Weights
#'
#' \code{kernel_weights} computes Gaussian kernel weights given a data matrix \code{X} and a scale parameter \code{phi}. Namely,
#' the lth weight \code{w[l]} is given by
#' \deqn{
#' w[l] = exp(-phi ||X[,i]-X[,j]||^2)
#' }, where the lth pair of nodes is (\code{i},\code{j}).
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param phi The nonnegative parameter that controls the scale of kernel weights
#' @useDynLib cvxbiclustr
#' @return A vector \cite{w} of weights for convex clustering.
kernel_weights <- function(X,phi=1) {
  storage.mode(X) <- "double"
  p <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  phi <- as.double(phi)
  w <- double(n*(n-1)/2)
  sol <- .C('kernel_weights',X=X,p=p,n=n,phi=phi,w=w)
  return(weights=sol$w)
}

#' Edge-Incidence Matrix of Weights Graph
#'
#' Construct the edge-incidence matrix of the weights graph.
#'
#' @param w Weights vector
#' @param n Number of points being clustered
#' @import Matrix
create_edge_incidence <- function(w,n) {
  P <- vec2tri(w@i+1,n)
  nEdges <- nrow(P)
  E <- Matrix(data=0,nrow=nEdges,ncol=n,sparse=TRUE)
  r <- 1:nEdges
  c <- P[,1]
  E[(c-1)*nEdges + r] <- 1
  c <- P[,2]
  E[(c-1)*nEdges + r] <- -1
  return(E)
}

#' Create adjacency matrix from V
#'
#' \code{create_adjacency} creates an n-by-n sparse adjacency matrix from the matrix of centroid differences.
#'
#' @param V Matrix of centroid differences
#' @param Phi Edge-incidence matrix
#' @import Matrix
#' @export
create_adjacency <- function(V,Phi) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0)
  n <- ncol(Phi)
  m <- length(connected_ix)
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)

  if (m > 0) {
    ix <- integer(m)
    jx <- integer(m)
    for (i in 1:m) {
      ix[i] <- which(Phi[connected_ix[i],]==1)
      jx[i] <- which(Phi[connected_ix[i],]==-1)
    }
    A[(jx-1)*n + ix] <- 1
  }
  return(A)
}

#' Find clusters
#'
#' \code{find_clusters} uses breadth-first search to identify the connected components of the corresponding
#' adjacency graph of the centroid differences vectors.
#'
#' @param A adjacency matrix
#' @export
#' @import igraph
find_clusters <- function(A) {
  G <- graph.adjacency(A, mode = 'upper')
  n <- nrow(A)
  node_seen <- logical(n)
  cluster <- integer(n)
  k <- 1
  for (i in 1:n) {
    if (!node_seen[i]) {
      connected_set <- graph.bfs(G, root=i, unreachable = FALSE)$order
      node_seen[connected_set] <- TRUE
      cluster[connected_set] <- k
      k <- k + 1
    }
  }
  nClusters <- k - 1
  size <- integer(nClusters)
  for (j in 1:nClusters) {
    size[j] <- length(which(cluster == j))
  }
  return(list(cluster=cluster, size=size))
}

#' Weights Graph Adjacency Matrix
#'
#' Constructs the adjacency matrix of the weights graph. This is useful to determine the connectivity of the weights graph.
#'
#' @param w Weights vector
#' @param n Number of points being clustered
#' @import Matrix
weights_graph <- function(w,n) {
  #  k <- which(w > 0)
  ix <- vec2tri(w@i+1,n)
  i <- ix[,1]
  j <- ix[,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}
