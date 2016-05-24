vec2Lfd <- function(bwtvec, rangeval=c(0,1))
{
#VEC2LFD converts a vector of length m to a linear differential
#  operator object of order m.  The range of the
#  functional data object in any cell is set to RANGEVAL.
#  In the event that BWTVEC is already a linear differential operator
#  object, it returns the object.

# Last modified 13 June 2013 by Spencer Graves
# Previously modified 10 December 2005

#  return BWTVEC if it is of class LFD

if (inherits(bwtvec, "Lfd")) {
    Lfdobj <- bwtvec
    return(Lfdobj)
}

#  check BWTVEC
if(!is.numeric(bwtvec))
    stop('first argument is not a linear differential operator ',
         'and is not numeric.')

isVec <- is.vector(bwtvec)
if(!isVec){
    bdim <- dim(bwtvec)
    if(is.null(bdim))
        stop('first argument is neither a linear differential operator ',
             'nor a vector (but is numeric)')
    if(length(bdim)>2)
        stop('first argument is neither a linear differential operator ',
             'nore a vector (but is an array of more than 2 dimensions).')
    if(!(1 %in% bdim))
        stop('first argument is neither a linear differential operator ',
             'nor a vector (but is a numeric matrix)')
    bwtvec <- as.numeric(bwtvec)
}

#if (!(is.vector(bwtvec) ||
#    stop("Argument not a vector and not a linear differential operator.")

m <- length(bwtvec)

#  set up the list object for the homogeneous part

if (m==0) {
    #  if derivative is zero, BWTLIST is empty
    bwtlist <- NULL
} else {
    conbasis <- create.constant.basis(rangeval)
    bwtlist  <- vector("list", m)
    for (j in 1:m) bwtlist[[j]] <- fd(bwtvec[j], conbasis)
}

#  define the Lfd object

Lfdobj <- Lfd(m, bwtlist)

return(Lfdobj)

}
