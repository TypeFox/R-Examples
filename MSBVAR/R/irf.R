"irf" <- function(varobj, nsteps, A0=NULL)
{
    if(inherits(varobj, "VAR")){
        return(irf.VAR(varobj, nsteps, A0=chol(varobj$mean.S)))
    }
    if(inherits(varobj, "BVAR")){
        return(irf.BVAR(varobj, nsteps, A0=chol(varobj$mean.S)))
    }
    if(inherits(varobj, "BSVAR")){
        return(irf.BSVAR(varobj, nsteps, A0=solve(varobj$A0.mode)))
    }
}

## "irf.msbsvar" <- function(gibbs, msbsvar, nsteps){
##     .Call("irf.msbsvar.cpp", gibbs, msbsvar, as.integer(nsteps))
## }


# 2010-08-06 : corrected classing mechanisms

"irf.VAR" <- function(varobj, nsteps, A0=chol(varobj$mean.S)){
    output <- .Call("irf.var.cpp", as.double(varobj$ar.coefs),
                    as.integer(dim(varobj$ar.coefs)),
                    as.integer(nsteps),
                    A0)

    attr(output, "class") <- c("irf", "irf.VAR")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"irf.BVAR" <- function(varobj, nsteps, A0=chol(varobj$mean.S))
{
    output <- irf.VAR(varobj, nsteps, A0=A0)
    attr(output, "class") <- c("irf", "irf.BVAR")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"irf.BSVAR" <- function(varobj, nsteps, A0=solve(varobj$A0.mode))
{
    output <- irf.VAR(varobj, nsteps, A0=A0)
    attr(output, "class") <- c("irf", "irf.BSVAR")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

# Updated 2006-01-29 to transpose rows and columns so the responses
# are in the rows.

# updated 2006-02-24 to use correct ranges for plots

"plot.irf" <- function(x, varnames = attr(x, "eqnames"), ...)
{
    #plot.irf.VAR(x, varnames = NULL, ...)
    if(inherits(x, "irf.VAR"))
    {
	plot.irf.VAR(x, varnames = varnames, ...)
    }
    if(inherits(x, "irf.BVAR"))
    {
	plot.irf.BVAR(x, varnames = varnames, ...)
    }
    if(inherits(x, "irf.BSVAR"))
    {
	plot.irf.BSVAR(x, varnames = varnames, ...)
    }
}

"plot.irf.VAR" <- function (x, varnames = NULL, ...)
{
    impulses <- x$mhat
    m <- dim(impulses)[1]
    nsteps <- dim(impulses)[3]
    stacked.impulses <- t(matrix(impulses, m^2, nsteps))
    stacked.impulses <- matrix(stacked.impulses, nsteps*m, m)
    dim(impulses) <- c((m^2), nsteps)
    impulses <- ts(t(impulses), start = c(0, 1), frequency = 1)
    minmax <- apply(stacked.impulses, 2, range)
    print(minmax)

    j <- 1
    par(mfrow = c(m, m), mai = c(0.25, 0.25, 0.15, 0.25), omi = c(0.15,
        0.75, 1, 0.15))
    for (i in 1:m^2) {

        plot(impulses[, i], xlab = "", ylab = "", ylim = c(minmax[,ceiling(i/m)]))
        abline(h = 0)
        if (i <= m) {
            mtext(varnames[i], side = 3, line = 2)
        }
        if ((i - 1)%%m == 0) {
            mtext(varnames[j], side = 2, line = 3)
            j <- j + 1
        }
    }
    mtext("Shock to", side = 3, line = 3, outer = T)
    mtext("Response in", side = 2, line = 3, , outer = T)
}

"plot.irf.BVAR" <- function (x, varnames = NULL, ...)
{
    plot.irf.VAR(x, varnames, ...)
}

"plot.irf.BSVAR" <- function (x, varnames = NULL, ...)
{
    plot.irf.VAR(x, varnames, ...)
}


## "irf.VAR.DEPRICATED" <- function(varobj, nsteps, A0=chol(varobj$mean.S))
##   {
##     ar.coef <- varobj$ar.coef
##     m<-dim(ar.coef)[1]                   # Capture the number of variables
##      p<-dim(ar.coef)[3]                   # Capture the number of lags
##      B<-array(0,c(m,m,max(nsteps,p)))     # Make an array for the AR coefs

##      for(l in 1:(min(p,nsteps)))
##        { B[,,l]<-ar.coef[,,l] }

##      mhat<-matrix(0,ncol=m*m,nrow=nsteps)
##      dim(mhat)<-c(m,m,nsteps)             # Make an array to hold IRF
##      mhat[,,1]<-A0                        # Identification condition
##                                           # for IRF
##      for(i in 2:nsteps)                   # Compute the IRF
##        { for(j in 1:(i-1))
##            { mhat[,,i]<-mhat[,,i] + (mhat[,,(i-j)]%*%B[,,j]) }
##        }
##     output <- list(B=B,mhat=mhat)
##     attr(output, "class") <- c("irf.VAR")
##     return(output)
## }
