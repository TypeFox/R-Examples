##function r = geneRank(W,ex,d)
##%function r = geneRank(W,ex,d)
##% 
##% GeneRank is a modification of the PageRank algorithm.
##% input data is   W: connectivity  matrix (zero/one, symmetric with zero diag)
##%                 ex: vector of expression levels (non-negative)
##%                 d: parameter in algorithm
##%
##% output is   r: vector of rankings
##%
##% March 09/2004
##%
##% Reference: GeneRank: Using search engine technology for the analysis
##%            of microarray experiments,       
##%            by Julie L. Morrison, Rainer Breitling, 
##%            Desmond J. Higham and David R. Gilbert, 
##%            submitted for publication.
##
##ex = abs(ex);
##norm_ex = ex/max(ex);
##w = sparse(W);
##degrees = sum(W);
##ind = find(degrees == 0);
##degrees(ind) = 1;
##D1 = sparse(diag(1./degrees));
##A = eye(size(w)) - d*(w'*D1); % eye produces identity matrix
##b = (1-d)*norm_ex;
##r = A\b;
##
##endfunction

geneRank <- function(W,ex,d, max.degree=Inf){

  ex = abs(ex)

  ## normalize expression values
  norm_ex = ex/max(ex)

  ## try sparse Matrices in R => later
  ##w = sparse(W)
  dimW = dim(W)[1]
  if(dim(W)[2]!=dimW) stop("W must be a square matrix.")
  
  ## get the in-degree for every node
  ## from KEGG we get a directed graph
  ## thus, the column sums correspond to
  ## to the in-degree of a particular gene
  degrees = pmin(max.degree, pmax(1, colSums(W), na.rm=T))

  ## A = Identity Matrix with dimensions
  ## same as the adjacency matrix
  A=matrix(0, nrow = dimW, ncol = dimW)
  diag(A) = 1
  
  ## produce a matrix with the degrees on
  ## the diagonal
  D1=matrix(0, nrow = dimW, ncol = dimW)
  diag(D1) = 1.0/degrees

  ## divide the in-degrees of the gene
  ## by the overall in-degree(colSum) of the gene
  ## => kind of normalizing the in-degrees
  A = A - d*(t(W) %*% D1)

  ## here, we give 1-d 'for free'
  b  = (1-d) * norm_ex

  ## we want to solve:
  ## (I - d W^t D^-1)r = (1-d)ex which is the Jacobi of the PageRank
  ## where A = (I - d W^t D^-1)
  ## and   b = (1-d)ex
  ## therefore:  Ar = b
  r = as.numeric(solve(A,b))
  return(r)
}
