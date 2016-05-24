#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Appell's F1 hypergeometric function
## 
## Time-stamp: <[appellf1.R] by DSB Don 19/04/2012 16:58 (CEST)>
##
## Description:
## Main function of the package.
##
## History:
## 26/03/2012   file creation
## 27/03/2012   - include debug option
##              - first finished version
## 28/03/2012   add reference for rkf45
## 30/03/2012   add reference for chyp.f
## 19/04/2012   add option to choose algorithm for hyp2f1 function
#####################################################################################



##' Compute Appell's F1 hypergeometric function
##'
##' This function is a wrapper for Fortran code written by F. D. Colavecchia and
##' G. Gasaneo, which is available at
##' \url{http://cpc.cs.qub.ac.uk/summaries/ADSJ}.  The corresponding background
##' reference is F. D. Colavecchia, G. Gasaneo and J. E. Miraglia (2001):
##' Numerical evaluation of Appell's F1 hypergeometric function, Computer
##' Physics Communications 138:29-43.
##'
##' External code in \dQuote{rkf45.f90} by L. F. Shampine and H. A. Watts is
##' used which is available from netlib.org at
##' \url{http://www.netlib.org/ode/rkf45.f}.  It is published in G. E. Forsythe,
##' M. A. Malcolm and C. B. Moler (1977): Computer Methods for Mathematical
##' Computations, Prentice-Hall.  Its performance is illustrated in F. Shampine,
##' H. A. Watts and S. Davenport (1976): Solving non-stiff ordinary differential
##' equations - the state of the art, SIAM Review 18:376-411.
##' 
##' The expert user can specify the actual computation with the parameter
##' \code{userflag}. Here the values 1 and 2 correspond to ODE integration and
##' series summation, respectively, while 0 decides between these two methods
##' based on the parameter values. Other possible values are 15-17 and 21-30,
##' each referring to an equation in F. D. Colavecchia et al. (2001). The
##' default value of \code{userflag} is -1 and leaves the algorithm decision to
##' the Fortran program, the result of which is returned in the list element
##' \code{algoflag}. Here the additional values 5 and 6 correspond to simple and
##' polynomial transformations, respectively.
##'
##' @param a complex parameter of Appell's F1  
##' @param b1 complex parameter of Appell's F1  
##' @param b2 complex parameter of Appell's F1  
##' @param c complex parameter of Appell's F1  
##' @param x numeric variable
##' @param y numeric variable
##' @param debug debug mode? (default is \code{FALSE})
##' @param userflag user flag variable (not used by default, expert option)
##' @param hyp2f1 which algorithm should be used for computing the Gaussian
##' hypergeometric function? See \code{\link{hyp2f1}} for details.
##'
##' @example examples/appellf1.R
##' 
##' @return A list with the algorithm and user flags as well as the 
##' complex value of Appell's F1 hypergeometric function.
##'
##' @export
##' @keywords math
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
appellf1 <-
    function(a,
             b1,
             b2,
             c,
             x,
             y,
             debug=FALSE,
             userflag=-1L,
             hyp2f1=c("michel.stoitsov", "forrey"))
{
    ## make sure all parameters are complex
    a <- as.complex(a)
    b1 <- as.complex(b1)
    b2 <- as.complex(b2)
    c <- as.complex(c)

    ## and variables are double
    x <- as.double(x)
    y <- as.double(y)

    ## debug must be logical
    debug <- as.logical(debug)

    ## make sure userflag is integer
    userflag <- as.integer(userflag)

    ## check value if it is present
    if(! (identical(userflag, -1L)))
    {       
        allowedFlags <- as.integer(c(0, # Depending on variables, either ODE
                                        # integr. or B&Ch Series
                                     1, # ODE Integr.
                                     2, # B&Ch Series
                                     15, 16, 17, 21, 22, 23, # Transformations
                                     24, 25, 26, 27, 28, 29, 30))
        if(! (userflag %in% allowedFlags))
            stop("user-specified flag not allowed!")
        else
            if(debug) cat("Using user-specified flag", userflag,
                          "for the computation\n") 
    }        
    
    ## allocate flag variable
    algoflag <- integer(1L)

    ## determine choice of algorithm for the Gaussian hypergeometric function
    hyp2f1 <- match.arg(hyp2f1)
    hyp2f1 <- switch(hyp2f1,
                     forrey=1L,
                     michel.stoitsov=2L)
    
    ## allocate return variable
    val <- complex(1L)
    
    ## then call Fortran to do the rest:
    results <- .Fortran(f1,
                        a=a,
                        b1=b1,
                        b2=b2,
                        c=c,
                        x=x,
                        y=y,
                        algoflag=algoflag,
                        userflag=userflag,
                        debug=debug,
                        val=val,
                        hyp2f1=hyp2f1)

    ## finally return everything
    return(results[c("algoflag", "userflag", "val")])
}
