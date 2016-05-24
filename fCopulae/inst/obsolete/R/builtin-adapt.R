

# Title: adapt -- multidimensional numerical integration
# Package: adapt
# Version: 1.0-4
# Author: FORTRAN by Alan Genz, 
#     S by Mike Meyer, R by Thomas Lumley and Martin Maechler
# Description: Adaptive Quadrature in up to 20 dimensions
# Depends:
# License: Unclear (Fortran) -- code in Statlib's ./S/adapt
# Maintainer: Thomas Lumley <tlumley@u.washington.edu>
# Packaged: Fri Apr 20 11:38:07 2007; thomas


# [from Statlib's original  http://lib.stat.cmu.edu/S/adapt ]
# This code contains an S function and supporting C and Fortran code for
# adaptive quadrature.  The underlyling fortran code is purported to
# work in from 2 to 20 dimensions.  The code is set up to dynamically
# load from a central library area.  If you can not do dynamic loading,
# you may need to build a staticly loaded version.  The adapt S function
# calls load.if.needed to do the dynamic loading.  You will have to
# change the functions used here (probably to call library.dynam).
# S code written by Michael Meyer (mikem@andrew.cmu.edu).
# October, 1989.


# 2002-03-14  Martin Maechler  <maechler@stat.math.ethz.ch>
# * DESCRIPTION (Version): 1.0-3 --> CRAN
# * R/adapt.R (adapt): use defaults for minpts, maxpts, eps;
#   more logical maxpts default (for ndim >= 7) using rulcls
# * man/adapt.Rd: extended example
# 2002-03-13  Martin Maechler  <maechler@stat.math.ethz.ch>
# * DESCRIPTION (Version): 1.0-2
# * man/adapt.Rd: indentation, using \code{.}, etc;
#   example also tries p=5 dimensions
# * R/adapt.R: clean up (spaces)
# 2002-01-09  Martin Maechler  <maechler@stat.math.ethz.ch>
# * R/adapt.R: do not use .Alias anymore
# 2001-06-29  Thomas Lumley <tlumley@u.washington.edu>
# * move (improved!) integrate() into base, using .Call() etc.


# Message-ID: <4AD7A74B.3020108@math.wsu.edu>
# Date: Thu, 15 Oct 2009 15:50:51 -0700
# From: Alan Genz <genz@math.wsu.edu>
# User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.21) 
#     Gecko/20090402 SeaMonkey/1.1.16
# MIME-Version: 1.0
# To: Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
# CC: Alan C Genz <alangenz@wsu.edu>
# Subject: Re: adapt
# References: <4AD3032B.4090801@itp.phys.ethz.ch>
# In-Reply-To: <4AD3032B.4090801@itp.phys.ethz.ch>
# Content-Type: text/plain; charset=ISO-8859-1; format=flowed
# Content-Transfer-Encoding: 7bit
# Status:  O
# Dear Prof. Wuertz,
# Thank you for your message and your interest in my adaptive integration
# Fortran code. I would be pleased if you included my code in your open
# source R fCopulae package under the Gnu GPL2 license. You have my
# permission to do this.
# Sincerely,
# Alan Genz


################################################################################


adapt <- function (ndim, lower, upper, minpts = 100, maxpts = NULL,
    functn, eps = 0.01, ...)
{
    keep.trying <- is.null(maxpts)

    if (ndim == 1) { ## fudge for 1-d functions
    warning("Using integrate() from base package for 1-d integration")
        if (keep.trying) maxpts <- minpts
    return(integrate(functn,lower,upper,subdivisions=maxpts,rel.tol=eps,...))
    }
    ## else ndim >= 2 :

    ## Check to make sure that upper and lower are reasonable lengths
    ## Both the upper and lower limits should be at least of length ndim
    if (length(lower) < ndim || length(upper) < ndim)#MM: dropped 'at least':
    stop(paste("The lower and upper vectors need to have ndim elements\n",
           "Your parameters are:  ndim", ndim, ", length(lower)",
           length(lower), ", length(upper)", length(upper), "\n"))
    ff <-
    if(length(list(...)) && length(formals(functn)) > 1)
        function(x) functn(x, ...)
    else functn # .Alias
    rulcls <- 2^ndim + 2*ndim^2 + 6*ndim + 1 #-> ../src/adapt.f

    ## maxpts should be large enough.  Prefer 10*rulclc, but use 2*rulclc.
    if (keep.trying)
        maxpts <- max(minpts, 500, 2 * rulcls)
    else {
        if (minpts >= maxpts) {
            warning(paste("maxpts must be > minpts.\n",
                          "Maxpts has be increased to  minpts + 1"))
            maxpts <- minpts + 1
        }
        ##
        if (maxpts < 2 * rulcls) {
            warning(paste(
                "You have maxpts (= ", maxpts, ") too small\n",
                "It needs to be at least 2 times 2^ndim + 2*ndim^2 + 6*ndim+1\n",
                "It has been reset to ", 2 * rulcls, "\n", sep=""))
            maxpts <- 2 * rulcls
        }
    }

    repeat {
    lenwrk <- (2*ndim + 3)* (1 + maxpts/rulcls)/2# mandated in adapt source

    x <- .C("cadapt",
        as.integer(ndim),
        as.double(lower),
        as.double(upper),
        minpts = as.integer(minpts),
        maxpts = as.integer(maxpts),
        ## now pass ff and current environment
        ff, rho = environment(),
        as.double(eps),
        relerr = double(1),
        lenwrk = as.integer(lenwrk),
        value = double(1),    # will contain the value of the integral
        ifail = integer(1),
        PACKAGE = "fCopulae")[
            c("value", "relerr", "minpts", "lenwrk", "ifail")]

    if (x$ifail == 1 && keep.trying)
        maxpts <- maxpts*2
    else
        break
    }
    if(x$ifail)
    warning(x$warn <-
        c("Ifail=1, maxpts was too small. Check the returned relerr!",
          paste("Ifail=2, lenwrk was too small. -- fix adapt() !\n",
            "Check the returned relerr!"),
          "Ifail=3: ndim > 20 -- rewrite the fortran code ;-) !",
          "Ifail=4, minpts > maxpts; should not happen!",
          "Ifail=5, internal non-convergence; should not happen!"
          )[x$ifail])

    class(x) <- "integration"
    
    x
}


# ------------------------------------------------------------------------------


print.integration <- function(x, ...) {
    print(noquote(sapply(x, format, ...)),...)
    invisible(x)
}


################################################################################

