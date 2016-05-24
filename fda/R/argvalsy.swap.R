argvalsy.swap = function(argvals=NULL, y=NULL, basisobj=NULL)
{
##
## 1.  if(is.null(y)) ...
##
  if(is.null(y)){
    if(is.null(argvals))stop("'y' is missing with no default")
#   Store argvals as y
    cat("'y' is missing, using 'argvals'\n") 
    y <- argvals
    argvals <- NULL 
  }
##
## 2.  if(is.null(argvals))argvals <- seq(basisobj$rangeval, dim(y)[1])
##
  dimy <- dim(as.array(y))
  if(is.null(argvals)){
    {
      if(is.null(basisobj))
        basisobj <- create.bspline.basis(basisobj)
      else {
        if(is.numeric(basisobj)) {
          if(length(basisobj)>1)
            basisobj <- create.bspline.basis(basisobj)
          else 
              basisobj <- create.bspline.basis(norder=basisobj)
        }
        else {
          if(inherits(basisobj, 'fd'))basisobj <- basisobj$basis
          else 
            if(inherits(basisobj, 'fdPar'))
              basisobj <- basisobj$fd$basis
        }
      }
    }
    a01 <- basisobj$rangeval
    if(is.null(a01))
      stop('basisobj does not have a required ',
           'rangeval component.')
#    
    n <- dimy[1]
    cat(paste("'argvals' is missing;  using seq(", a01[1],
            ", ", a01[2], ", length=", n, ")\n"))       
    argvals <- seq(a01[1], a01[2], length=n)
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
##
## 3.  if(length(dim(argvals)) == length(dim(y))) ... 
##
  dima <- dim(as.array(argvals))
  {
    if(length(dimy) == length(dima)){
      if(any(dimy != dima))
        stop("dimensions of 'argvals' and 'y' must be compatible;\n",
             "  dim(argvals) = ", paste(dima, collapse=' x '),
             ";  dim(y) = ", paste(dimy, collapse=' x ') )
#     Check basisobj
      {
        if(inherits(basisobj, 'fd'))basisobj <- basisobj$basis
        else {
          if(inherits(basisobj, 'fdPar'))
            basisobj <- basisobj$fd$basis
          else {
            if(inherits(basisobj, 'array')){
              fd. <- fd(basisobj)
              basisobj <- fd.$basis
            }
            else { 
              if(inherits(basisobj, 'integer'))
                basisobj <- create.bspline.basis(argvals, norder=basisobj)
              else {
                if(is.null(basisobj))
                  basisobj <- create.bspline.basis(argvals)
                else
                  if(!inherits(basisobj, 'basisfd'))
                    stop("'basisobj' is NOT a functional basis",
                         " object (class 'basisfd');  class = ",
                         class(basisobj)[1])
              }
            }
          }
        }
      }
      a01 <- basisobj$rangeval
      arng <- range(argvals)
      if((a01[1]<=arng[1]) && (arng[2]<=a01[2]))
        return(list(argvals=argvals, y=y, basisobj=basisobj))
#
      yrng <- range(y)
      if((a01[1]<=yrng[1]) && (yrng[2]<=a01[2])){
        cat(paste("'argvals' is NOT contained in basisobj$rangeval",
                ", but 'y' is;  swapping 'argvals' and 'y'.\n"))
        return(list(argvals=y, y=argvals, basisobj=basisobj)) 
      }
#      
      stop("Neither 'argvals' nor 'y' are contained in ",
           "basisobj$rangeval")
    }
  }        
##
## 4.  If(length(dimy) < length(dima)) swap ...
##
  if(length(dimy)<length(dima)){
    cat(paste("Swapping 'y' and 'argvals', because 'y' is ",
            "simpler,\n  and 'argvals' should be;  now ",
            "dim(argvals) = ", paste(dimy, collapse=" x "),
            ";  dim(y) = ", paste(dima, collapse=" x "),"\n" )) 
    y. <- argvals
    argvals <- y
    y <- y.
#
    d. <- dima
    dima <- dimy
    dimy <- d.
  }   
#
  if(any(dima != dimy[1:length(dima)]))
    stop("A dimension of 'argvals' does not match 'y':\n",
         "  dim(argvals) = ", paste(dima, collapse=" x "),
         ";  dim(y) = ", paste(dimy, collapse=" x ") )      
##        
## 5.  Check compatibility of argvals with basisobj
##        
  {
    if(inherits(basisobj, 'fd'))basisobj <- basisobj$basis
    else {
      if(inherits(basisobj, 'fdPar'))
        basisobj <- basisobj$fd$basis
      else {
        if(inherits(basisobj, 'array')){
          fd. <- fd(basisobj)
          basisobj <- fd.$basis
        }
        else { 
          if(inherits(basisobj, 'integer'))
            basisobj <- create.bspline.basis(argvals, norder=basisobj)
          else {
            if(is.null(basisobj))
              basisobj <- create.bspline.basis(argvals)
            else
              if(!inherits(basisobj, 'basisfd'))
                stop("'basisobj' is NOT a functional basis",
                     " object (class 'basisfd');  class = ",
                     class(basisobj)[1])
          }
        }
      }
    }
  }
  a01 <- basisobj$rangeval
  arng <- range(argvals)
  if((a01[1]<=arng[1]) && (arng[2]<=a01[2]))
    return(list(argvals=argvals, y=y, basisobj=basisobj))
#
  stop("'argvals' are not contained in basisobj$rangeval")
}
