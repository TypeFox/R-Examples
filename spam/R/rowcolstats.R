# This is file ../spam/R/rowcolstats.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



rowSums.spam <- function(x,...) {
  return( .Fortran("rowsums",
                   as.double(x@entries), as.integer(x@colindices), as.integer(x@rowpointers),
                   x@dimension[1],
                   rs=vector("double",x@dimension[1]),
                   NAOK=.Spam$NAOK, PACKAGE="spam")$rs)
  
}

colSums.spam <- function(x,...) {
  return( .Fortran("colsums",
                   as.double(x@entries), as.integer(x@colindices), as.integer(x@rowpointers),
                   x@dimension[1],
                   cs=vector("double",x@dimension[2]),
                   NAOK=.Spam$NAOK, PACKAGE="spam")$cs)
}

rowMeans.spam <- function(x,...) {
  return( .Fortran("rowmeans",
                   as.double(x@entries), as.integer(x@colindices), as.integer(x@rowpointers),
                   x@dimension[1],x@dimension[2],
                   as.logical(.Spam$structurebased),
                   rm=vector("double",x@dimension[1]),
                   NAOK=.Spam$NAOK, PACKAGE="spam")$rm)
}

colMeans.spam <- function(x,...) {
   return( .Fortran("colmeans",
                as.double(x@entries), as.integer(x@colindices), as.integer(x@rowpointers),
                x@dimension[1],x@dimension[2],
                as.logical(.Spam$structurebased),
                cm=vector("double",x@dimension[2]),vector("integer",x@dimension[2]),
                NAOK=.Spam$NAOK, PACKAGE="spam")$cm)
}



setMethod("rowSums","spam",rowSums.spam)
setMethod("colSums","spam",colSums.spam)
setMethod("rowMeans","spam",rowMeans.spam)
setMethod("colMeans","spam",colMeans.spam)
