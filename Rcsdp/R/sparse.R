## A simple class for sparse symmetric (triplet) matrices.
##
## Mostly copied from sparse.R in Kurt Hornik's relations package

simple_triplet_sym_matrix <-
function(i, j, v, n = max(c(i,j)),check.ind=FALSE)
{
  if (check.ind & any(i<j)) {
      stop("Index arguments 'i' and 'j' do not point to the lower triangle. Swapping indices")
    }
  structure(list(i = i, j = j, v = v, n = n),
            class = "simple_triplet_sym_matrix")
}

as.simple_triplet_sym_matrix <-
function(x,...)
    UseMethod("as.simple_triplet_sym_matrix")

as.simple_triplet_sym_matrix.simple_triplet_sym_matrix <- function(x,...) x

as.simple_triplet_sym_matrix.matrix <-
function(x,check.sym=FALSE)
{
    if(prod(dim(x)) == 0L)
        return(simple_triplet_sym_matrix(integer(), integer(), c(x), n=nrow(x)))

    if (nrow(x) != ncol(x))
      stop("Argument 'x' must be a square matrix")

    if (check.sym && (sum(abs(x-t(x))/2) != 0))
        stop("Argument 'x' must be a symmetric matrix")
    
    ind <- which(x != vector(typeof(x), 1L), arr.ind = TRUE)
    if (length(ind) == 0)
      return(.simple_triplet_zero_sym_matrix(nrow(x)))
             
    ind <- ind[ind[,1L] >= ind[,2L],,drop=FALSE]
    simple_triplet_sym_matrix(ind[, 1L], ind[, 2L], x[ind],
                          n = nrow(x))
}

as.matrix.simple_triplet_sym_matrix <-
function(x, ...)
{
    n <- x$n
    y <- matrix(vector(typeof(x$v), n * n), n, n)
    ind <- cbind(x$i, x$j)
    y[ind] <- x$v
    offdiag <- x$i != x$j
    y[ind[offdiag,c(2L,1L),drop=FALSE]] <- x$v[offdiag]
    y
}

as.vector.simple_triplet_sym_matrix <-
  function(x, ...)
  {
    as.vector(as.matrix(x,...))
  }

as.double.simple_triplet_sym_matrix <-
  function(x,...)
  {
    as.double(as.vector(x,...))
  }

dim.simple_triplet_sym_matrix <-
function(x)
    c(x$n, x$n)

all.equal.simple_triplet_sym_matrix <-
  function(target,current,...)
  {
    if(!inherits(current,"simple_triplet_sym_matrix"))
      stop("Argument 'current' is not of class 'simple_triplet_sym_matrix'")
    
    if(!inherits(target,"simple_triplet_sym_matrix"))
       stop("Argument 'target' is not of class 'simple_triplet_sym_matrix'")

    cur.ord <- order(current$j,current$i)
    current$i <- current$i[cur.ord]
    current$j <- current$j[cur.ord]
    current$v <- current$v[cur.ord]
    
    targ.ord <- order(target$j,target$i)
    target$i <- target$i[targ.ord]
    target$j <- target$j[targ.ord]
    target$v <- target$v[targ.ord]

    all.equal.list(current,target,...)
  }

## Utitilies for creating special simple triplet matrices:

.simple_triplet_zero_sym_matrix <-
function(n, mode = "double")
    simple_triplet_sym_matrix(integer(), integer(), vector(mode, 0L),n)

.simple_triplet_diag_sym_matrix <-
function(x, n)
{
    x <- rep(x, length.out = n)
    i <- seq_len(n)
    simple_triplet_sym_matrix(i, i, x, n)
}

.simple_triplet_random_sym_matrix <-
function(n,occ=.1,nnz=occ*n*(n+1)/2,rfun=rnorm,seed=NULL,...)
  {
    if(!missing(seed) & !is.null(seed))
      set.seed(seed)

    # sample indices in the lower triangle
    iind <- sample(n*(n+1)/2,nnz)
    ind <- 1:(n*(n+1)/2) + sapply(rep(1:n,times=n:1), function(x) x*(x-1)/2)
    ind <- ind[iind]
    simple_triplet_sym_matrix(i=((ind-1) %% n) + 1,
                          j=((ind-1) %/% n) + 1,
                          v=rfun(nnz,...),
                          n=n)
  }

