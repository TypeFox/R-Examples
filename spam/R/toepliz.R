# This is file ../spam/R/toepliz.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     





########################################################################

"circulant.spam" <- function(x, n=NULL, eps = .Spam$eps)
{
  if (!(is.vector(x)|is.list(x)) )
    stop("'x' is not a vector or a list")
  
  if( is.list(x)) {
    if (!identical(length(x),2L))
      stop("Argument 'x' needs to be a list with two elements")
    if (is.null(n))
      stop("'n' needs to be given")
    ind <- x[[1]]
    x <- x[[2]]
    sel <- (ind <= n)&(abs(x)>eps)
    ind <- ind[sel]
    x <- x[sel]
   
  }else{
    n <- length(x)
    ind <- (1:n)[abs(x) > eps]
    x <- x[ind]
  }
  
  n <- as.integer(n)
  len <- as.integer(length( ind)[1]) # see ?length@value
  if(identical(len,0L))
    return(new("spam", rowpointers = c(1L, rep.int(2L, n)), 
               dimension = as.integer(c(n, n))))
#      subroutine circulant(nrow,len, x,j, a,ja,ia)
  nz <- n*len
  z <- .Fortran('circulant',n, len,
                as.double(x), as.integer(ind),
                entries= vector("double", nz),
                colindices = vector("integer", nz),
                rowpointers = vector("integer",  n + 1),
                 NAOK = .Spam$NAOK, PACKAGE = "spam")

                
  newx <- new("spam")
  slot(newx, "entries", check = FALSE) <- z$entries
  slot(newx, "colindices", check = FALSE) <- z$colindices
  slot(newx, "rowpointers", check = FALSE) <- z$rowpointers
  slot(newx, "dimension", check = FALSE) <- c(n, n)
  return(newx)

}

toeplitz.spam <- function(x,y=NULL, eps = .Spam$eps)
{
  if (!is.vector(x)) 
    stop("'x' is not a vector")
  n <- length(x)

  if (!is.null(y)){
    if (!identical(length(y),n))
      stop("Length of 'y' and 'x' do not match")
    fullx <- c(rev(y[-1]),x)
  } else {
    fullx <- c(rev(x[-1]),x)
  }

  ind <- (1:(2*n-1))[abs(fullx) > eps]
  fullx <- fullx[ind]

 
  n <- as.integer(n)
  len <- as.integer(length( ind)[1]) # see ?length@value
  if(identical(len,0L))
    return(new("spam", rowpointers = c(1L, rep.int(2L, n)), 
               dimension = as.integer(c(n, n))))
#      subroutine toeplitz(nrow,len, x,j, a,ja,ia,kk)
  nz <- n*len
  z <- .Fortran('toeplitz',n, len,
                as.double(fullx), as.integer(ind),
                entries= vector("double", nz),
                colindices = vector("integer", nz),
                rowpointers = vector("integer",  n + 1),
                nnz=as.integer(1),
                NAOK = .Spam$NAOK, PACKAGE = "spam")

                
  newx <- new("spam")
  slot(newx, "entries", check = FALSE) <- z$entries[1:z$nnz]
  slot(newx, "colindices", check = FALSE) <- z$colindices[1:z$nnz]
  slot(newx, "rowpointers", check = FALSE) <- z$rowpointers
  slot(newx, "dimension", check = FALSE) <- c(n, n)
  return(newx)


}


#
