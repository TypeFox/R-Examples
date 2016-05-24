# Function from package 'fda' (c) 2014

int2Lfd <- function(m=0)
{
#INT2LFD converts a nonnegative integer to a linear differential
#  operator object that is equivalent to D^m.  The range of the
#  functional data object in any cell is set to [0,1], and is
#  not actually used when a linear differential operator object
#  of this nature is applied.  
#  In the event that m is already a linear differential operator
#  object, it returns the object immediately.  Thus, INT2LFD can
#  be used to screen whether an object is an integer or not.

#  Last modified 17 September 2005

#  check M

if (inherits(m, "Lfd")) {
    Lfdobj <- m
    return(Lfdobj)
}

if (!is.numeric(m)) 
    stop("Argument not numeric and not a linear differential operator.")
 

if (length(m) != 1) stop("Argument is not a scalar.")

if (round(m) != m)  stop("Argument is not an integer.")

if (m < 0)   stop("Argument is negative.")

#  all the checks passed, set up a functional data object
#  The range will not be used in this case, and can be set
#  to [0, 1]

#  set up the list object for the homogeneous part

if (m==0) {
    #  if derivative is zero, BWTLIST is empty
    bwtlist <- NULL
} else {
    basisobj <- create.constant.basis(c(0,1))
    bwtlist  <- vector("list", m)
    for (j in 1:m) bwtlist[[j]] <- fd(0, basisobj)
}

#  define the Lfd object

Lfdobj <- Lfd(m, bwtlist)

return(Lfdobj)

}
