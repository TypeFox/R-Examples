# functions to support sparse matrices from the Matrix package
.not.dsTMatrix <- function(x,cl)
  {
    if(!inherits(x,cl))
      stop("Argument 'x' must be of class ", cl)
    as.simple_triplet_sym_matrix(as(x,"dsTMatrix"))
  }

as.simple_triplet_sym_matrix.dgeMatrix <- function(x) .not.dsTMatrix(x,"dgeMatrix")
as.simple_triplet_sym_matrix.dsyMatrix <- function(x) .not.dsTMatrix(x,"dsyMatrix")
as.simple_triplet_sym_matrix.dpoMatrix <- function(x) .not.dsTMatrix(x,"dpoMatrix")
as.simple_triplet_sym_matrix.dgTMatrix <- function(x) .not.dsTMatrix(x,"dgTMatrix")
as.simple_triplet_sym_matrix.dsCMatrix <- function(x) .not.dsTMatrix(x,"dsCMatrix")
as.simple_triplet_sym_matrix.dsTMatrix <-
  function(x)
  {
    if(!inherits(x,"dsTMatrix"))
      stop("Argument 'x' must be of class 'dsTMatrix'")

    # simple_triplet_sym_matrix store entries in the
    # lower triangle only
    if (x@uplo %in% c("u","U")) {
      x@uplo <- "l"
      tmp <- x@i
      x@i <- x@j
      x@j <- tmp
    }
    
    x <- as(x,"dsTMatrix")
    simple_triplet_sym_matrix(i=x@i+1,
                              j=x@j+1,
                              v=x@x,
                              n=nrow(x))
  }

#setAs("simple_triplet_sym_matrix","dsTMatrix",
#  function(x)
#  {
#    new("dgTMatrix",i=as.integer(x$i-1),j=as.integer(x$j-1),x=v,Dim=as.integer(c(x$n,x$n)),uplo="L")
#  })

as.simple_triplet_sym_matrix.ddiMatrix <- function(x)
  {
    if(!inherits(x,"ddiMatrix"))
       stop("Argument 'x' must be of class 'ddiMatrix'")

    v <- if (x@diag=="U") rep(1,nrow(x)) else x@x
    simple_triplet_sym_matrix(i=1:nrow(x),
                              j=1:nrow(x),
                              v=v,
                              n=nrow(x))
  }

