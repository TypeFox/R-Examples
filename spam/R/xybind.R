# This is file ../spam/R/xybind.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     









########################################################################
"rbind.spam" <-
function(...,deparse.level=0)
{
  if (deparse.level!=0) warning("Only 'deparse.level=0' implemented, coerced to zero,")
  addnargs <- ifelse(missing(deparse.level),0,1)

  nargs <- nargs()-addnargs
  if (nargs == 0)     return( NULL)
  args <- list(...)
  if (!is.null( names( args)))  {
    warning("Names of arguments are ignored")
    names( args) <- NULL
  }
  args[which(sapply(args, is.null))] <- NULL

  # nargs needs an update
  nargs <- length(args) -  addnargs
  

  if (nargs == 0)     return( NULL)
  if (nargs == 1)     return( args[[1]])
  if (nargs == 2) {
    # we distinguish between the cases:
    #    1 spam, spam
    #    2 spam, numeric (scalar, vector, matrix)
    #    3 numeric, spam
    #    4 numeric, numeric
    
    # Case 1: this is the quick way
    if( is.spam(args[[1]]) & is.spam(args[[2]])) {

      if(ncol(args[[1]])!=ncol(args[[2]]))
        stop("Arguments have differing numbers of columns, in rbind.spam()",call.=FALSE)
    
      nrow1 <- args[[1]]@dimension[1]

      newx <- new("spam")
      newx@entries <- c(args[[1]]@entries, args[[2]]@entries)
      newx@colindices <- c(args[[1]]@colindices,  args[[2]]@colindices)
      newx@rowpointers <- c(args[[1]]@rowpointers,
                            args[[2]]@rowpointers[-1]+args[[1]]@rowpointers[nrow1+1]-as.integer(1))
      newx@dimension <- c(nrow1+args[[2]]@dimension[1],args[[1]]@dimension[2])
      return(newx)
    }
    # Case 2:  spam, numeric (scalar, vector, matrix)
    #    if scalar, coherce it first to vector of appropriate length,
    #    if vector, attach dimension.
    if( is.spam(args[[1]]) & is.numeric(args[[2]])) {
      Xdim <- args[[1]]@dimension
      Ylen <- length(args[[2]])
      if (Ylen==1) {
        Xlen <- Xdim[2]
        args[[2]] <- rep( args[[2]], Xlen)
        dim( args[[2]]) <- c(1,Ylen)

      } else   if (is.vector(args[[2]]))
        dim(args[[2]]) <- if( Xdim[1]==1) c(Ylen,1) else c(1,Ylen)
      Ydim <- dim(args[[2]])

      if(Xdim[2]!=Ydim[2])
        stop("Arguments have differing numbers of columns, in rbind.spam()",call.=FALSE)
    
 
      newx <- new("spam")
      newx@entries <- c(args[[1]]@entries, as.double(t(args[[2]])))
      newx@colindices <- c(args[[1]]@colindices,  rep.int(as.integer(1:Ydim[2]),Ydim[1]))
      newx@rowpointers <- c(args[[1]]@rowpointers,
                            seq.int(args[[1]]@rowpointers[Xdim[1]+1], by=Ydim[2], length.out=Ydim[1]+1)[-1])
      newx@dimension <- c(Xdim[1]+Ydim[1],Ydim[2])
      return(newx)
    }
    # Case 3:  numeric (scalar, vector, matrix), spam
    #    similar as above
    if( is.numeric(args[[1]]) & is.spam(args[[2]])) {
      Xlen <- length( args[[1]])
      Ydim <- args[[2]]@dimension
      if (Xlen==1) {
        Xlen <- Ydim[2]
        args[[1]] <- rep( args[[1]], Xlen)
        dim( args[[1]]) <- c(1,Xlen)
      } else   if (is.vector(args[[1]]))
        dim(args[[1]]) <- if ( Ydim[1]==1) c(Xlen,1) else  c(1,Xlen)
      Xdim <- dim(args[[1]])

      if(ncol(args[[2]])!=Xdim[2])
        stop("Arguments have differing numbers of columns, in rbind.spam()",call.=FALSE)
    

      newx <- new("spam")
      newx@entries <- c(as.double(t(args[[1]])), args[[2]]@entries )
      newx@colindices <- c(rep.int(as.integer(1:Xdim[2]),Xdim[1]),
                           args[[2]]@colindices)
      newx@rowpointers <- c(seq.int(1, by=Xdim[2], length.out=Xdim[1]),
                            args[[2]]@rowpointers + Xlen)
      newx@dimension <- c(Ydim[1]+Xdim[1],Ydim[2])
      return(newx)
    }
    # Case 4: numeric,numeric
    #    result is a cleaned spam object.
    if( is.numeric(args[[1]]) & is.numeric(args[[2]]))
      return( as.spam.matrix( rbind(args[[1]],args[[2]]))  )

    stop("Not all argument are of class 'spam' and 'numeric', in rbind.spam()",
         call.=FALSE)
  } else {
    # "recursive" approach only, e.g. no checking
    tmp <- rbind.spam( args[[1]],args[[2]])
    for ( i in 3:nargs)
      tmp <- rbind.spam( tmp,args[[i]])
    return( tmp)
  }
}
  
  
"cbind.spam" <-
function(...,deparse.level=0)
{
  if (deparse.level!=0) warning("Only 'deparse.level=0' implemented, coerced to zero,")
  addnargs <- ifelse(missing(deparse.level),0,1)

  nargs <- nargs()-addnargs
  if (nargs == 0)     return( NULL)
  args <- list(...)
  if (!is.null( names( args)))  {
    warning("Names of arguments are ignored")
    names( args) <- NULL
  }
  args[which(sapply(args, is.null))] <- NULL
  nargs <- length(args) -  addnargs
  if (nargs == 0)     return( NULL)
  if (nargs == 1)     return( args[[1]])
  if (nargs == 2) {

    Ydim <- if (is.spam(args[[2]])) args[[2]]@dimension  else
          dim(args[[2]])

    if (is.numeric(args[[1]])) {
      # we do _not_ have a spam object 
      if (is.vector(args[[1]]))  {
        if (is.null(Ydim)) {
          # if Ydim is NULL, then Y is a vector as well and we need special treatment.
          args[[1]] <- spam.numeric( args[[1]],
                                    max( length(args[[1]]),length(args[[2]])))      
                                    
        } else {
          # "standard" case, a vector (scalar) with a matrix
          args[[1]] <- spam.numeric( args[[1]], Ydim[1])                  # scalar
        }
      } else {
        # we have a regular matrix
        args[[1]] <- as.spam.matrix(args[[1]])
      }
    } else if (!is.spam(args[[1]])) {
      # we have anything
      stop("Not all argument are of class 'spam' and 'numeric', in cbind.spam()",
           call.=FALSE)
    }
      
    Xdim <- args[[1]]@dimension

    # now, X=args[[1]] is a spam, the treatment of Y is easier:
    if (is.numeric(args[[2]])) {
      # we do _not_ have a spam object 
      if (is.vector(args[[2]]))    args[[2]] <- spam.numeric( args[[2]], Xdim[1])  
      else  args[[2]] <- as.spam.matrix(args[[2]])
    } else if (!is.spam(args[[2]])) {
      stop("Not all argument are of class 'spam' and 'numeric', in cbind.spam()",
           call.=FALSE)
    }
    
    Ydim <- args[[2]]@dimension

    if(Xdim[1]!=Ydim[1])
         stop("Arguments have differing numbers of rows, in cbind.spam()",call.=FALSE)

    XYlen <- args[[1]]@rowpointers[Xdim[1]+1]+args[[2]]@rowpointers[Xdim[1]+1]-2L
    z <- .Fortran("cbind", Xdim[2], Xdim[1], Ydim[2], XYlen,
                  args[[1]]@entries, args[[1]]@colindices, args[[1]]@rowpointers,
                  args[[2]]@entries, args[[2]]@colindices, args[[2]]@rowpointers,
                  entries=vector( "double", XYlen),
                  colindices=vector( "integer", XYlen),
                  rowpointers=vector( "integer", Xdim[1]+1),
                  NAOK=.Spam$NAOK,PACKAGE = "spam")

    if (FALSE) {
    # a loop would be (in R...):
      for (i in 1:nrow) {
        if (args[[1]]@rowpointers[i]<args[[1]]@rowpointers[i+1])
          stend1 <- args[[1]]@rowpointers[i]:(args[[1]]@rowpointers[i+1]-1)
        else stend1 <- NULL
        if (args[[2]]@rowpointers[i]<args[[2]]@rowpointers[i+1])
          stend2 <- args[[2]]@rowpointers[i]:(args[[2]]@rowpointers[i+1]-1)
        else stend2 <- NULL
        entries <- c( entries, args[[1]]@entries[stend1], args[[2]]@entries[stend2])
        colindices <-  c(  colindices, args[[1]]@colindices[stend1], args[[2]]@colindices[stend2]+Xdim[2])
      }
    }
    newx <- new("spam")
    slot(newx,"entries", check=FALSE)     <- z$entries
    slot(newx,"colindices", check=FALSE)  <- z$colindices
    slot(newx,"rowpointers", check=FALSE) <- z$rowpointers
    slot(newx,"dimension", check=FALSE)   <- c(Ydim[1],Xdim[2]+Ydim[2])
    return(newx)
    
  } else {
    # "recursive" approach only, e.g. no checking
    tmp <- cbind.spam( args[[1]],args[[2]])
    for ( i in 3:nargs)
      tmp <- cbind.spam( tmp,args[[i]])
    return( tmp)
  }
}
  
  


setMethod("rbind","spam",rbind.spam)
setMethod("cbind","spam",cbind.spam)


