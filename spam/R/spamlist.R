# This is file ../spam/R/spamlist.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


"spam.list" <-  function(x, nrow = 1, ncol = 1,eps = .Spam$eps) {
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  if (!is.list(x)|(length(x)<2)|(length(x)>3))
    stop("Argument 'x' needs to be a list with two or three elements")
  # two cases: list of length
  # -  two (matrix with two columns called ind* and the elements)
  # -  three (each one column called i*, j*.   

  if (identical(length(x),2L)) {
    indnr <- pmatch("ind",names(x)) 
    if (is.na(indnr)) stop("Argument 'x' needs an element called 'indices'")
    elenr <- ifelse( identical( indnr,1L), 2L, 1L)
    
    nz <- length( x[[elenr]])

    dimx <- dim(x[[indnr]])
    if (is.null(dimx)||(dimx[2] != 2))  stop("Indices should have two columns")
    if (dimx[1] != nz) stop("Number of indices does not match with number of elements")
    
    ir <- as.integer(x[[indnr]][,1])
    jc <- as.integer(x[[indnr]][,2])
  } else {
    inr <- pmatch("i",names(x)) 
    jnr <- pmatch("j",names(x))
    if (is.na(inr)||is.na(jnr)) stop("Argument 'x' needs elements called 'i' and 'j'")
    elenr <- c(1:3)[-c(inr,jnr)]
    nz <- length( x[[elenr]])
    
    ir <- as.integer(x[[inr]])
    jc <- as.integer(x[[jnr]])

    if ((length(ir) != nz)||(length(jc) != nz))
        stop("Number of indices does not match with number of elements")
  }
  if (identical(nz, 0L))
    return(new("spam",rowpointers=c(1L,rep.int(2L,as.integer(nrow))), dimension=as.integer(c(nrow,ncol))))
  if (any( ir <= 0) || any( jc <= 0))
      stop("Indices need to be positive")
  if (any(!is.finite(x[[elenr]]))) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[[elenr]][!is.finite(x[[elenr]])] <- 0
  }
  nrow <- as.integer(ifelse(missing(nrow),max(ir),nrow))
  ncol <- as.integer(ifelse(missing(ncol),max(jc),ncol))
  z <- .Fortran(ifelse(toupper(.Spam$listmethod=='PE'),"triplet3csr","triplet2csr"),
                nrow=nrow, ncol=ncol,
                nz=nz,
                as.double(x[[elenr]]),ir,jc,
                entries=vector("double",nz),
                colindices=vector("integer",nz),
                rowpointers=vector("integer",nrow+1),eps,
                NAOK=TRUE, PACKAGE = "spam"
                )
#  print(z)
  if (identical(z$nz, 0L))
    return(new("spam",rowpointers=c(1L,rep.int(2L,nrow)), dimension=c(nrow,ncol)))
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$entries[1:z$nz]
  slot(newx,"colindices",check=FALSE) <- z$colindices[1:z$nz]
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newx,"dimension",check=FALSE) <- c(nrow,ncol)
  
  return(newx)
}

setMethod("as.spam","list",   { function(x,eps) spam.list(x,eps=eps)})
setMethod("spam","list",spam.list)
