#' BLACS Sums
#' 
#' Sum across a process grid.
#' 
#' For advanced users only.
#' 
#' @param ICTXT
#' BLACS ICTXT.
#' @param SCOPE
#' Rows, cols, or both.
#' @param m,n
#' Problem size.
#' @param x
#' Local values.
#' @param lda
#' Leading dimension.
#' @param RDEST
#' Row destination.
#' @param CDEST
#' Col destination.
#' 
#' @rdname blacs-sums
#' @export
base.igsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call(R_igsum2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

#' @rdname blacs-sums
#' @export
base.dgsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgsum2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}



#' BLACS Max
#' 
#' Max value across a process grid.
#' 
#' For advanced users only.
#' 
#' @param ICTXT
#' BLACS ICTXT.
#' @param SCOPE
#' Rows, cols, or both.
#' @param m,n
#' Problem size.
#' @param x
#' Local values.
#' @param lda
#' Leading dimension.
#' @param RDEST
#' Row destination.
#' @param CDEST
#' Col destination.
#' 
#' @rdname blacs-max
#' @export
base.igamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call(R_igamx2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

#' @rdname blacs-min
#' @export
base.dgamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgamx2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}



#' BLACS Min
#' 
#' Min value across a process grid.
#' 
#' For advanced users only.
#' 
#' @param ICTXT
#' BLACS ICTXT.
#' @param SCOPE
#' Rows, cols, or both.
#' @param m,n
#' Problem size.
#' @param x
#' Local values.
#' @param lda
#' Leading dimension.
#' @param RDEST
#' Row destination.
#' @param CDEST
#' Col destination.
#' 
#' @rdname blacs-min
#' @export
base.igamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call(R_igamn2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

#' @rdname blacs-min
#' @export
base.dgamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgamn2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}


#' BLACS Point to Poin
#' 
#' Sent value across a process grid.
#' 
#' For advanced users only.
#' 
#' @param ICTXT
#' BLACS ICTXT.
#' @param SCOPE
#' Rows, cols, or both.
#' @param m,n
#' Problem size.
#' @param x
#' Local values.
#' @param lda
#' Leading dimension.
#' @param RDEST
#' Row destination.
#' @param CDEST
#' Col destination.
#' 
#' @rdname blacs-p2p
#' @export
base.dgesd2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgesd2d1, as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

#' @rdname blacs-p2p
#' @export
base.dgerv2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    pbdMPI::comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgerv2d1, as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}


