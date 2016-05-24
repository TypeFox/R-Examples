# This is file ../spam/R/tcrossprod.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

crossprod.spam <- function(x, y=NULL) {
    dimx <- dim(x)
    if( is.null(y)) {
        if(!is.spam(x)) return(crossprod(x))

        if(dimx[2]==1L) return(matrix( sum(x@entries^2)))

        return( .spam.matmul(t.spam(x),x))
    }
    if( (!is.spam(x)) & (!is.spam(y))) return(crossprod(x,y))
    
    return( .spam.matmul(t(x),y))
    
}
tcrossprod.spam <- function(x, y=NULL) {
    dimx <- dim(x)
    if( is.null(y)) {
        if(!is.spam(x)) return(tcrossprod(x))

        if(dimx[2]==1L) return(matrix( sum(x@entries^2)))

        return( .spam.matmul(x,t.spam(x)))
    }
    if( (!is.spam(x)) & (!is.spam(y))) return(tcrossprod(x,y))
    
    return( .spam.matmul(x,t(y)))
    
}


setMethod("crossprod",signature(x="spam",y="missing"), crossprod.spam)
setMethod("crossprod",signature(x="ANY",y="spam"), crossprod.spam)
setMethod("tcrossprod",signature(x="spam",y="missing"), tcrossprod.spam)
setMethod("tcrossprod",signature(x="ANY",y="spam"), tcrossprod.spam)
