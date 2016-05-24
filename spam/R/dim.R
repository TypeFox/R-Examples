# This is file ../spam/R/dim.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

# This is the actual dim...

"dim<-.spam" <- function(x, value) {
    if (is.spam(x)) {
        
        dimx <- x@dimension
        pdim <- prod(dimx)
        vlen <- prod(value)
        if( !identical(pdim,vlen))
            stop( sprintf("dims [product %d] do not match the length of object [%d]. Do you want `pad`",
                          pdim,vlen))
        
        if (length(value)>2)
            stop("dims should be of length 1 or 2")
        if (identical(length(value),1L))
            return( c(x) )

        if(any(dimx<1))
            stop("the dims contain negative values")
        
        tmp <- cbind(st=rep(1:dim(x)[1],diff(x@rowpointers)), nd=x@colindices)
        ind <- tmp[,1]+(tmp[,2]-1)*dimx[1] - 1

        slist <- list(i = ind%%value[1]   +1,
                       j = ind%/%value[1] +1,
                       x@entries)
        
        return( spam.list( slist, nrow=value[1], ncol=value[2],
                          eps = .Machine$double.eps))
        
        
    } else  {
        dim(x) <- value
        x
    }
}


########################################################################
# dim and derivatives

"pad<-.spam" <- function(x,value) {
  if ( (min(value)<1 ) || any(!is.finite(value)))
    stop("dims should be postive integers.")
  if (!identical( length(value), 2L)) stop("dims should be of length 2.")
  dimx <- x@dimension
  last <- value[1]+1

  # In three steps:
  #  1) Address col truncation
            # to safe time, we also take into account if we have fewer or equal rows
  #  2) Augment rows
  #  3) if fewer rows and more columns, truncate
  # In any case, dimensions are fixed at the end.
  
  # If fewer cols required, we run reducedim
  if (dimx[2]>value[2]){
#     subroutine reducedim(a,ja,ia,eps,bnrow,bncol,k,b,jb,ib)
    z <- .Fortran("reducedim",
                  oldra=as.double(x@entries),
                  oldja=x@colindices,
                  oldia=x@rowpointers,
                  eps=.Spam$eps,
                  as.integer(min(value[1],dimx[1])),as.integer(value[2]),
                  nz=1L,
                  entries=vector("double",length(x@entries)),
                  colindices=vector("integer",length(x@entries)),
                  rowpointers=vector("integer",last),
                  NAOK = .Spam$NAOK, PACKAGE = "spam")
    if (identical(z$nz,1L) )
      return(new("spam",rowpointers=c(1L,rep.int(2L,as.integer(value[1]))),
                 dimension=as.integer(value)))

    nz <- z$nz-1
    slot(x,"entries",check=FALSE) <- z$entries[1:nz]
    slot(x,"colindices",check=FALSE) <- z$colindices[1:nz]
    slot(x,"rowpointers",check=FALSE) <- z$rowpointers[1:min(last,dimx[1]+1)]
  }
  # augment rows
  if  (dimx[1]<value[1]){
    slot(x,"rowpointers",check=FALSE) <- c(x@rowpointers,rep.int(
                           x@rowpointers[length(x@rowpointers)],value[1]-dimx[1]))
  }
  # special case: fewer rows and more columns, truncate
  if ((dimx[1]>value[1])&(dimx[2]<value[2])) {
    lastelement <- (x@rowpointers[last]-1)
    slot(x,"entries",check=FALSE) <- x@entries[1:lastelement]
    slot(x,"colindices",check=FALSE) <- x@colindices[1:lastelement]
    slot(x,"rowpointers",check=FALSE) <- x@rowpointers[1:last]
  }
        
  slot(x,"dimension",check=FALSE) <- as.integer(value)                  
  return(x)

}



setMethod("dim",   "spam", function(x) x@dimension )
setMethod("dim<-",   "spam", get("dim<-.spam"))

setGeneric("pad<-", function(x, value) standardGeneric("pad<-"))
setMethod("pad<-",   "spam", get("pad<-.spam"))
setMethod("pad<-",   "matrix",
          function(x, value) {
              if (!identical( length(value), 2L)) stop("dims should be of length 2.")
              tmp <- matrix(0, value)
              mr <- 1:min(value[1], nrow(x))
              mc <- 1:min(value[2], ncol(x))
              tmp[mr,mc] <- x[mr,mc]
              return(tmp)
          })
