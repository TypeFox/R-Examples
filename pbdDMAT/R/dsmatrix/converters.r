### TODO put in appropriate as.___.r


dmat_ldim <- function(nrows, rank=comm.rank()) # FIXME add communicator
{
  rem <- nrows %% comm.size()
  
  n <- as.integer(nrows / comm.size())
  
  if (rank < rem)
    n <- n + 1L
  
  return( n )
}

# starting index
dmat_index <- function(nrows)
{
  if (comm.rank() == 0)
    start <- 1L
  else
  {
    cs <- comm.size() - 2L
    chunks <- sapply(0L:cs, dmat_ldim, nrows=nrows)
    start <- sum(chunks[1L:comm.rank()]) + 1L
  }
  
  return( start )
}



# -----------------------------------------------------------
# x = dsmatrix
# -----------------------------------------------------------

setMethod("as.dmat", signature(x="dsmatrix"),
  function(x)
  {
    Data <- convert_csr_to_dense(dim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    llb <- new("dmat", Data=Data, dim=x@dim, ldim=x@ldim, storage="llb")
    
    return( llb )
  }
)



setMethod("as.dsvector", signature(x="dsmatrix"),
  function(x)
  {
    if (x@dim[2L] != 1)
      comm.stop("not yet supported")
    
    y <- new("dsvector", length=x@dim[1L], llength=x@ldim[1L], Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind, storage=x@storage)
    
    return( y )
  }
)



# -----------------------------------------------------------
# x = dmat
# -----------------------------------------------------------

setMethod("as.dsmatrix", signature(x="matrix"),
  function(x)
    as.dsmatrix(as.dmat(x))
)



setMethod("as.dsmatrix", signature(x="dmat"),
  function(x)
  {
    l <- convert_dense_to_csr(x@Data)
    sparse <- new("dsmatrix", Data=l$Data, dim=x@dim, ldim=x@ldim, row_ptr=l$row_ptr, col_ind=l$col_ind, storage="csr")
    
    return( sparse )
  }
)




# -----------------------------------------------------------
# x = matrix
# -----------------------------------------------------------



setMethod("as.dmat", signature(x="matrix"),
  function(x)
  {
    dim <- dim(x)
    nrows <- dim[1L]
    
    nrows.local <- dmat_ldim(nrows)
    ldim <- c(nrows.local, dim[2L])
    
    start <- dmat_index(nrows)
    end <- start + nrows.local - 1L
    
    if (nrows.local == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- x[start:end, ]
    
    dmat <- new("dmat", Data=Data, dim=dim, ldim=ldim)
    
    return( dmat )
  }
)


