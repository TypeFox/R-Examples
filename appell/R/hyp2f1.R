#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Appell's F1 hypergeometric function
## 
## Time-stamp: <[hyp2f1.R] by DSB Don 19/04/2012 17:00 (CEST)>
##
## Description:
## Compute the Gaussian hypergeometric function with complex arguments.
##
## History:
## 18/04/2012   file creation
#####################################################################################



##' Compute the Gaussian hypergeometric function with complex arguments
##'
##' Two different algorithms can be used.
##' 
##' The first, default, algorithm uses Fortran code in \dQuote{hyp_2F1.f90} from
##' N. L. J.  Michel and M. V. Stoitsov, which is available at
##' \url{http://cpc.cs.qub.ac.uk/summaries/AEAE}. The corresponding background
##' reference is N. L. J. Michel and M. V. Stoitsov (2008): Fast computation of
##' the Gauss hypergeometric function with all its parameters complex with
##' application to the Pöschl-Teller-Ginocchio potential wave functions,
##' Computer Physics Communications 178:535-551. 
##' 
##' The second algorithm uses Fortran code in \dQuote{cyp.f} from R. C. Forrey
##' is used which is available at \url{http://physics.bk.psu.edu/codes/chyp.f}.
##' The corresponding background reference is R. C. Forrey (1997): Computing the
##' hypergeometric function, Journal of Computational Physics 137:79-100.
##' 
##' @param a complex parameter
##' @param b complex parameter
##' @param c complex parameter
##' @param z complex variable
##' @param algorithm either \dQuote{michel.stoitsov} (default) or
##' \dQuote{forrey} (see the details)
##'
##' @example examples/hyp2f1.R
##' 
##' @return The complex value of the Gaussian hypergeometric function.
##'
##' @export
##' @keywords math
##' @encoding utf-8
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
hyp2f1 <-
    function(a,
             b,
             c,
             z,
             algorithm=c("michel.stoitsov", "forrey"))
{
    ## make sure all parameters are complex
    a <- as.complex(a)
    b <- as.complex(b)
    c <- as.complex(c)
    z <- as.complex(z)

    ## determine choice of algorithm for the Gaussian hypergeometric function
    algorithm <- match.arg(algorithm)
    algorithm <- switch(algorithm,
                        forrey=1L,
                        michel.stoitsov=2L)
    
    ## allocate return variable
    val <- complex(1L)
    
    ## then call Fortran to do the rest:
    results <- .Fortran(f21_sub,
                        a=a,
                        b=b,
                        c=c,
                        z=z,
                        hyp2f1=algorithm,
                        val=val)

    ## finally return the result
    return(results$val)
}
