setMethod("t", signature(x="dsmatrix"),
  function(x)
  {
    m <- x@dim[2L]
    n <- x@dim[1L]
    
    new_dim <- c(m, n)
    new_ldim <- c(dmat_ldim(m), n)
    
    out <- petsc_mattranspose(dim=x@dim, ldim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    
    if (new_ldim[1L] == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- out$Data
    
    ret <- new("dsmatrix", Data=Data, dim=new_dim, ldim=new_ldim, row_ptr=out$row_ptr, col_ind=out$col_ind)
    
    return( ret )
  }
)



setMethod("%*%", signature(x="dsmatrix", y="dsmatrix"),
  function(x, y)
  {
    # Matrix multiplication dimension crap
    if (x@dim[2L] != y@dim[1L])
      comm.stop("Error : non-conformable arguments.")
    
    m <- x@dim[1L]
    n <- y@dim[2L]
    
    new_dim <- c(m, n)
    new_ldim <- c(dmat_ldim(m), n)
    
    
    # Compute product
    out <- petsc_matmatmult(A_dim=x@dim, A_ldim=x@ldim, A_data=x@Data, A_row_ptr=x@row_ptr, A_col_ind=x@col_ind, 
                            B_dim=y@dim, B_ldim=y@ldim, B_data=y@Data, B_row_ptr=y@row_ptr, B_col_ind=y@col_ind)
    
    
    # wrangle return
    if (new_ldim[1L] == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- out$Data
    
    ret <- new("dsmatrix", dim=new_dim, ldim=new_ldim, Data=Data, row_ptr=out$row_ptr, col_ind=out$col_ind)
    
    return( ret )
  }
)



setMethod("solve", signature(a="dsmatrix"),
  function(a)
  {
    out <- petsc_invert(dim=a@dim, ldim=a@ldim, Data=a@Data, row_ptr=a@row_ptr, col_ind=a@col_ind)
    
    if (ldim[1L] == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- out$Data
    
    ret <- new("dmat", dim=a@dim, ldim=a@ldim, Data=Data)
    
    return( ret )
  }
)


