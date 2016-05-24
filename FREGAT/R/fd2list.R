# Function from package 'fda' (c) 2014

fd2list <- function(fdobj)
{
#  FD2LIST converts a univariate functional data object to a list
#  object, mainly for purposes of defining a linear differential
#  operator object.

#  For example, this code sets up a harmonic acceleration Lfd object
#    over the interval [0,365] for the daily weather data.
#  Lbasis  = create.constant.basis(c(0,365));  #  create a constant basis
#  Lcoef   = matrix(c(0,(2*pi/365)^2,0),1,3)   #  set up three coefficients
#  wfdobj  = fd(Lcoef,Lbasis)      # define an FD object for weight functions
#  wfdlist = fd2list(wfdobj)       # convert the FD object to a cell object
#  harmaccelLfd = Lfd(3, wfdlist)  #  define the operator object

#  Last modified 26 October 2005

#  get the coefficient matrix and the basis

    coef     <- fdobj$coefs
    coefsize <- dim(coef)
    nrep     <- coefsize[2]

#  check whether FDOBJ is univariate

    if (length(coefsize) > 2)
    	stop("FDOBJ is not univariate.")

    fdlist <- vector("list",0)
    for (i in 1:nrep) fdlist[[i]] <- fdobj[i]
    return(fdlist)
}
