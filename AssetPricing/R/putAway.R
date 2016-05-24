putAway <- function(odeRslt,type,jmax,qmax,soltype,x=NULL,prices=NULL) {
#
# Equation solved; now do a bunch of housekeeping to store
# the results in a convenient manner.  Note than odeRslt is
# a matrix of dimension
#            (a) nout x (1 + np + 2*qmax) --- for continuous pricing,
# and a smooth price sensitivity function, or
#            (b) nout x (1 + qmax + np + qmax) --- for discrete pricing or
# a piecewise linear price sensitivity function, or
#            (c) nout x (1 + 2*qmax) --- for a given pricing policy
# where we are solving only for the expected value v.
# In the foregoing np is equal to the total number of price functions.
# This is given by np = jmax*(qmax - (jmax-1)/2)
#
# The first column of odeRslt consists of the times at which the
# solution to the differential equation is evaluated.  In case (a)
# columns 2 through np+1 contain the corresponding values of the
# optimal price x, and the last qmax columns contain the corresponding
# expected values.  In case (b), columns 2 through qmax+1 contain
# the optimal expected values, and the last np columns contain the
# corresponding optimal prices.  In case (c) the last qmax columns
# contain expected values (corresponding to the prices in the given
# pricing policy).
#
# If type == "sip" then when putAway() is called, jmax is effectively
# 1. In the process of solving the differential equation, jmax is
# the length of the vector of group size probabilities.  If type is
# "dip" then there is a price for each (q,j) combination for j = 1,
# ..., min{jmax,q} and hence jmax has an impact on the number of
# price functions. However when type is "sip" then there is only
# one price for each q value and so the number of price functions
# is qmax, the same as it would be if jmax were equal to 1.

if(type=="sip") jmax <- 1
tvec <- odeRslt[,1]
tmax <- max(tvec)
np   <- jmax*(qmax - (jmax-1)/2)
nc   <- ncol(odeRslt)

# Put away the prices.
#
# If putAway() is being called by vsolve() then these prices are
# *not* necessarily optimal. In this case "x" is passed as an argument
# and nothing needs to be done to it.

# Note that if type == "dip" then the entry x[[i]] is equal to
# the function x_qj(t) where i = (j-1)*(qmax - j/2) + q.
if(is.null(x)) {
    if(soltype=="cont")
        xx <- as.data.frame(odeRslt[,2:(np+1)])
    else
        xx <- as.data.frame(odeRslt[,(2+qmax):(1+qmax+np)])
    if(soltype=="disc") {
        x <- lapply(xx,function(x,t){approxfun(x=t,y=x,method="constant",
                                    yleft=NA,yright=NA)},t=tvec)
        x <- lapply(x,a2sf,tt=tvec)
    } else x <- lapply(xx,function(x,t){splinefun(x=t,y=x)},t=tvec)

    ylim <- range(xx)
    attr(x,'tlim') <- c(0,tmax)
    attr(x,'ylim') <- ylim
    attr(x,'qmax') <- qmax
    attr(x,'jmax') <- jmax
    if(!is.null(prices)) attr(x,'prices') <- prices
    comment(x) <- "Optimal prices."
    class(x) <- "flap"
    if(soltype=="disc")
        class(x) <- c(class(x),"pwc.flap")  # pwc <--> piecewise constant.
    if(type=="dip") class(x) <- c(class(x),"di.flap") # di  <--> doubly indexed.
}

# Put away the optimal expected values.
if(soltype=="cont")
    vv <- as.data.frame(odeRslt[,(2+np):(1+np+qmax)])
else
    vv <- as.data.frame(odeRslt[,2:(1+qmax)])
v    <- lapply(vv,function(x,t){splinefun(x=t,y=x)},t=tvec)
ylim <- range(vv)
attr(v,'tlim') <- c(0,tmax)
attr(v,'ylim') <- ylim
attr(v,'qmax') <- qmax
attr(v,'jmax') <- 1
class(v) <- "flap"

# Put away the derivatives of the optimal expected values.
vd <- as.data.frame(odeRslt[,(nc-qmax+1):nc])
vdot <- lapply(vd,function(x,t){splinefun(x=t,y=x)},t=tvec)
ylim <- range(vd)
attr(vdot,'tlim') <- c(0,tmax)
attr(vdot,'ylim') <- ylim
attr(vdot,'qmax') <- qmax
attr(vdot,'jmax') <- 1
class(vdot) <- "flap"
rslt <- list(x=x,v=v,vdot=vdot)
class(rslt) <- "AssetPricing"
rslt
}
