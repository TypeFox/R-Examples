# This is file ../spam/R/kronecker.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     









kronecker.spam <- function(X,Y,FUN = "*", make.dimnames = FALSE, ...)
{
  if (make.dimnames) 
    warning("dimnames not supported within sparse matrices")
  lenx <- length(X)
  leny <- length(Y)
  kronlen <- lenx*leny
  
  if(is.spam(X)){    
    Xdim <- X@dimension
    Xentries <- X@entries
    Xcol <- X@colindices
    Xrow <- X@rowpointers
  }else if (is.vector(X)){
    Xentries <- X
    Xcol <- rep.int(as.integer(1),lenx)
    Xdim <- as.integer(c(lenx,1))
    Xrow <- seq_len(lenx+1)
  } else {
    Xentries <- as.double(t(X))
    Xdim <- dim(X)
    Xcol <- rep.int(as.integer(1:Xdim[2]),Xdim[1])
    Xrow <- seq.int(1, by=Xdim[2], length.out=Xdim[1]+1)
  }

  if(is.spam(Y)){
    Ydim <- Y@dimension
    Yentries <- Y@entries
    Ycol <- Y@colindices
    Yrow <- Y@rowpointers
  } else if (is.vector(Y)){
    Yentries <- Y
    Ycol <- rep.int(as.integer(1),leny)
    Ydim <- as.integer(c(leny,1))
    Yrow <- seq_len(leny+1)
  } else {
    Yentries <- as.double(t(Y))
    Ydim <- dim(Y)
    Ycol <- rep.int(as.integer(1:Ydim[2]),Ydim[1])
    Yrow <- seq.int(1, by=Ydim[2], length.out=Ydim[1]+1)
  }
  kronxy <- new("spam")

  if (FUN=="*") {
    z <- .Fortran("kroneckermult",Xdim[1],Xentries,Xcol,Xrow,
                  Ydim[1],Ydim[2],Yentries,Ycol,Yrow,
                  entries=vector( "double", kronlen),
                  colindices=vector( "integer", kronlen),
                  rowpointers=vector( "integer",Xdim[1]*Ydim[1]+1),
                  NAOK=.Spam$NAOK,PACKAGE = "spam")
    slot(kronxy, "entries", check=FALSE) <-  z$entries
   }else {
    z <- .Fortran("kronecker",Xdim[1],Xentries,Xcol,Xrow,
                  Ydim[1],Ydim[2],Yentries,Ycol,Yrow,
                  ent1=vector( "double",kronlen),ent2=vector( "double",kronlen),
                  colindices=vector( "integer",kronlen),
                  rowpointers=vector( "integer",Xdim[1]*Ydim[1]+1),
                  NAOK=.Spam$NAOK,PACKAGE = "spam")
    FUN <- match.fun(FUN)
    slot(kronxy, "entries", check=FALSE) <-  FUN(z$ent1,z$ent2,...)
    if (z$rowpointers[Xdim[1]*Ydim[1]+1]-1 < prod(Xdim,Ydim))
      warning("Sparseness structure of 'kronecker(X,Y)' preseved when applying 'FUN'.", call. = FALSE)
  }
  slot(kronxy, "colindices", check=FALSE) <- z$colindices
  slot(kronxy, "rowpointers", check=FALSE) <- z$rowpointers
  slot(kronxy, "dimension", check=FALSE) <- Xdim*Ydim   # no need to use prod

  return(kronxy)
}


"kronecker.default" <- base::kronecker
setGeneric("kronecker")

setMethod("kronecker","spam", kronecker.spam)
setMethod("kronecker",signature(X="spam",Y="spam"), kronecker.spam)
setMethod("kronecker",signature(X="spam",Y="ANY"), kronecker.spam)
setMethod("kronecker",signature(X="ANY",Y="spam"), kronecker.spam)
