# This is file ../spam/R/constructors.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

"rowpointers<-" <- function(x, value) {
    dimx <- x@dimension
    nnz1 <- x@rowpointers[dimx[1]+1]
    diffvalue <- diff(value)
    if ( any(!is.finite(value)))
        stop("row pointers should be postive integers.")
    if (!identical( length(x@rowpointers), length(value)))
        stop("wrong length of row pointers in `rowpointers<-`.", call.=FALSE)
    if (any(diffvalue<0))
        stop("row pointers are not monotone increasing in `rowpointers<-`.", call.=FALSE)
    if (any(diffvalue>dimx[2]))
        stop("row pointers have too large leaps in `rowpointers<-`.", call.=FALSE)
    if (value[1]<1)
        stop("first element of row pointers is < 1 in `rowpointers<-`.", call.=FALSE)
    if(value[dimx[1]+1] != nnz1)
        stop("last element of row pointers does not conform in `rowpointers<-`.", call.=FALSE)
       
    x@rowpointers <- as.integer(value)
    x
}

"colindices<-" <- function(x, value) {
    dimx <- x@dimension
    if ( any(!is.finite(value)))
        stop("column indices should be postive integers in `colindices<-`.", call.=FALSE)
    
    if ( any(value<1) | any(value> dimx[2]))
        stop("column indices exceed dimension `colindices<-`.", call.=FALSE)
    diffcolindices <- diff(value)     # positive values within each row
    if (all(diff(x@rowpointers)>1) && length(diffcolindices)>0)
                                        # only if we have multiple values
      if (identical( dimx[1], 1L)) {
        if   ( any(diffcolindices<1))
          stop("column indices are not ordered `colindices<-`.", call.=FALSE)
      } else {
        if ( any(diffcolindices[-(x@rowpointers[2:dimx[1]]-1)]<1))
          stop("column indices are not ordered `colindices<-`.", call.=FALSE)
      }
    x@colindices <- as.integer(value)
    x
}

"entries<-" <- function(x, value) {
    if (!identical( length(x@entries), length(value)))
        stop("wrong length in `entries<-`.", call.=FALSE)
    if (!.Spam$NAOK) {
        if (any(!is.finite(value)))
            stop("'NA/NaN/Inf' not allowed in `entries<-`.", call.=FALSE)
    }
    if (!is.numeric( value))
        stop("numerical required in `entries<-`.", call.=FALSE)
    
    x@entries <- as.double(value)
    x
}    

"dimension<-" <- function(x, value) {
    stop("modification through `dim' or `pad`", call.=FALSE)
}


