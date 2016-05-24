# This is file ../spam/R/permutation.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








checkpivot <- function(pivot, len, type="Pivot") {
  if(is.null(pivot))                       return()
  if(!is.vector(pivot))                        stop(paste(type,"is not a vector."))     
  pivot <- as.vector(pivot,"integer")
  if (!identical(length(pivot),len))           stop(paste(type,"of wrong length."))
  tmp <- sort.int(pivot)
  if(tmp[1]!=1 ||  any(tmp-seq_len(len)!=0))   stop(paste("Invalid",type))
  return()
}



"permutation.spam" <- function(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE){
  # eliminated .Internal calls as this creates a 'Note' on CRAN checks.
  # Only 1-2% timing loss, see end of the file.
  
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]

  if (is.null(P)&is.null(Q))
    stop("At least one permutation should be specified")

  nz <- A@rowpointers[nrow+1]-1


  if (check){
    checkpivot(P,nrow,"Permutation")
    checkpivot(Q,ncol,"Permutation")
  }
  
  if (is.null(Q)) {
#      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm)
#      B = P A
    P <- as.integer(P)
    if(ind) P <- order(P)
#    if(ind) P <- .Internal(order(T,F,P))
   
    z <- .Fortran("rperm",
                  nrow,
                  A@entries,A@colindices,A@rowpointers,
                  entries = vector("double",nz), 
                  colindices = vector("integer", nz),
                  rowpointers = vector("integer", nrow + 1),P, 
                  NAOK = .Spam$NAOK, PACKAGE = "spam")
  } else {  
    if (is.null(P)){
#      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,iwork) 
#      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), iwork(*)
#      double precision a(*), ao(*) 
#      B = A Q 
      Q <- as.integer(Q)
      if(ind) Q <- order(Q)
#       if(ind) Q <- .Internal(order(T,F,Q))
        z <- .Fortran("cperm",
                      nrow,
                      A@entries,A@colindices,A@rowpointers,
                      entries = vector("double",nz),
                      colindices = vector("integer", nz),
                      rowpointers = vector("integer", nrow + 1),
                      Q,
                      NAOK = .Spam$NAOK, PACKAGE = "spam")
      } else {  
#      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,pperm,qperm,iwork)
#      B = P A Q 
        Q <- as.integer(Q)
#        if(ind) Q <- .Internal(order(T,F,Q))
        if(ind) Q <- order(Q)
        P <- as.integer(P)
#        if(ind) P <- .Internal(order(T,F,P))
        if(ind) P <- order(P)
        z <- .Fortran("dperm",
                      nrow,
                      A@entries,A@colindices,A@rowpointers,
                      entries = vector("double",nz),
                      colindices = vector("integer", nz),
                      rowpointers = vector("integer", nrow + 1),
                      P,Q,
                      NAOK = .Spam$NAOK, PACKAGE = "spam")
      }   
    
  }
  newx <- new("spam")
  slot(newx, "entries", check = FALSE) <- z$entries
  slot(newx, "colindices", check = FALSE) <- z$colindices
  slot(newx, "rowpointers", check = FALSE) <- z$rowpointers
  slot(newx, "dimension", check = FALSE) <- c(nrow,ncol)
  return(newx)
}

permutation.matrix <- function(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE){
  nrow <- dim(A)[1]
  ncol <- dim(A)[1]
  
  if (is.null(P)&is.null(Q))     stop("At least one permutation should be specified")
  
  if (check){
    checkpivot(P,nrow,"Permutation")
    checkpivot(Q,ncol,"Permutation")
  }

  if (ind) {
    if (is.null(Q))     return(A[P,])
    if (is.null(P))     return(A[,Q])
    return(A[P,Q])
  } else {
    if (is.null(Q))     return(A[order(P),])
    if (is.null(P))     return(A[,order(Q)])
    return(A[order(P),order(Q)])
  }
}

setGeneric("permutation",function(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE)standardGeneric("permutation"))
setMethod("permutation","matrix",permutation.matrix)
setMethod("permutation","spam",permutation.spam)


### ss <- sample(1:100000)
### system.time( for( i in 1:1000) tt<-order(ss))
### system.time( for( i in 1:1000) tt<-.Internal(order(T,F,ss)))
