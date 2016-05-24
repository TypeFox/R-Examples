##' Solve for IRR (internal rate of return) or YTM (yield to maturity)
##' 
##' This function computes the internal rate of return at which the
##' net present value equals zero. It requires as input a function
##' that computes the net present value of a series of cash flows for
##' a given interest rate as well as the derivative of the npv with
##' respect to the interest rate (10,000 times this derivative is the
##' PVBP or DV01).  In this package, \code{irr.solve} is primarily
##' intended to be called by the \code{\link{irr}} and
##' \code{\link{bond.yield}} functions. It is made available for those
##' who want to find irr of more complex instruments.
##' 
##' The function \code{irr.solve} is basically an interface to the
##' general root finder \code{\link{newton.raphson.root}}. However, if
##' \code{\link{newton.raphson.root}} fails, \code{irr.solve} makes an
##' attempt to find the root using \code{\link{uniroot}} from the R
##' stats package. Uniroot uses bisection and it requires the root to
##' be bracketed (the function must be of opposite sign at the two end
##' points - lower and upper).
##' 
##'
##' @param f The function whose zero is to be found. An R function
##'     object that takes one numeric argument and returns a list of
##'     two components (value and gradient). In the IRR applications,
##'     these two components will be the NPV and its derivative
##' @param interval The interval c(lower, upper) within which to
##'     search for the IRR
##' @param r.guess The starting value (guess) from which the solver
##'     starts searching for the IRR
##' @param toler The argument \code{toler} to
##'     \code{\link{newton.raphson.root}}.  The IRR is regarded as
##'     correct if abs(NPV) is less than \code{toler}.  Otherwise the
##'     \code{irr.solve} returns \code{NA}
##' @param convergence The argument \code{convergence} to
##'     \code{\link{newton.raphson.root}}.
##' @param max.iter The maximum number of iterations of the
##'     Newton-Raphson procedure
##' @param method The root finding method to be used. The
##'     \code{default} is to try Newton-Raphson method
##'     (\code{\link{newton.raphson.root}}) and if that fails to try
##'     bisection (\code{\link{bisection.root}}). The other two
##'     choices (\code{newton} and \code{bisection} force only one of
##'     the methods to be tried.
##' @return The function irr.solve returns \code{NA} if the irr/ytm
##'     could not be found. Otherwise it returns the irr/ytm. When
##'     \code{NA} is returned, a warning message is printed
##' 
##' @author Prof. Jayanth R. Varma \email{jrvarma@@iimahd.ernet.in}
##' @export
irr.solve <- function(f, interval = NULL, r.guess = NULL, toler = 1e-6,
                      convergence = 1e-8, max.iter = 100,
                      method = c('default', 'newton', 'bisection')){
    ## Bisection needs an interval. If interval was not provided, we
    ## start with a wide interval starting with the theoretical bound
    ## of -1 (or -100 percent) on the lower side and the square root
    ## of the largest machine integer on the upper side. With 32 bit
    ## integers, the upper limit is over 40,000 or over 4 million
    ## percent. Even with 16 bit integers, the upper limit is over
    ## 18,000 percent.
    lower <- if (! is.null(interval)) max(interval[1], -1) else -1
    upper <- if (! is.null(interval) && is.finite(interval[2]))
                 interval[2] else sqrt(.Machine$integer.max)
    if(is.null(r.guess)){
        r.guess <- if(is.null(interval)) 0 else mean(interval)
    }
    if(lower > r.guess || upper < r.guess){
        warning('Guess value must be inside intervals')
        return(NA)
    }
    method <- match.arg(method)
    if(method != 'bisection'){
        ## Try Newton Raphson method
        res <- NA
        tryCatch(
            res <- newton.raphson.root(f=f, guess=r.guess, lower=lower,
                                       upper=upper, max.iter=max.iter,
                                       convergence=convergence,
                                       toler=toler),
            warning = function(war){}
        )
        if (! is.na(res)){
            return (res)
        }
    }
    if(method != 'newton'){
        ## Try bisection. 
        
        ## We add one to everything to get rid of possible negative
        ## lower limits and make everything positive as required by
        ## bisection.root to allow geometric steps. One is actually
        ## slightly more than 1.0 to make things strictly positive
        one <- 1.01
        ## Correspondingly we subtract one before calling the function
        fv <- function(x) f(x - one)$value
        res2 <- NA
        tryCatch(
            res2 <- bisection.root(f = fv, guess = r.guess + one,
                                   lower = lower + one,
                                   upper = upper + one, toler=toler),
            warning = function(war){}
        )
        if (! is.na(res2)){
            return (res2 - one)
        }
    }
    warning(.irr.warning.msg(method))
    return(NA)
}

## internal function not exported
.irr.warning.msg <- function(method){
    wmsg <- list(
        'default' = "Both Newton-Raphson and Bisection failed",
        'newton' = "Newton-Raphson failed",
        'bisection' = "Bisection failed")
    paste(wmsg[[method]], " to find IRR/YTM")
    ## cat('irr.warning with method = ', method, '\n')
    ## "IRR failed"
}
##' Find zero of a function by bracketing the zero and then using
##' bisection.
##'
##' Tries to find the zero of a function by using the bisection method
##' (\code{\link[stats]{uniroot}}). To call
##' \code{\link[stats]{uniroot}}, the zero must be bracketed by
##' finding two points at which the function value has opposite
##' signs. The main code in this function is a grid search to find
##' such a pair of points. A geometric grid of points between
##' \code{lower} and \code{guess} and also between \code{guess} and
##' \code{upper}. This grid is searched for two neighbouring points
##' across which the function changes sign. This brackets the root,
##' and then we try to locate the root by calling
##' \code{\link[stats]{uniroot}}
##' @param f The function whose zero is to be found. An R function
##'     object that takes one numeric argument and returns a numeric
##'     value. In an IRR application, this will be the NPV
##'     function. In an implied volatility application, the value will
##'     be the option price.
##' @param guess The starting value (guess) from which the solver
##'     starts searching for the root. Must be positive.
##' @param lower The lower end of the interval within which to search
##'     for the root. Must be positive.
##' @param upper The upper end of the interval within which to search
##'     for the root. Must be positive.
##' @param nstep THe number of steps in the grid search to bracket the
##'     zero. See details.
##' @param toler The criterion to determine whether a zero has been
##'     found. This is passed on to \code{\link[stats]{uniroot}}
##' @return The root (or NA if the method fails)
##' @export
##' @author Prof. Jayanth R. Varma
bisection.root <- function(f, guess, lower, upper, nstep = 100,
                           toler = 1e-6){
    if(lower <= 0 || upper <= 0 || guess <= 0){
        warning("lower, upper and guess must all be positive")
        return(NA)
    }
    if(lower > guess || upper < guess){
        warning('Guess value must be between lower and upper')
        return(NA)
    }
    ## We move from guess to lower in nsteps geometric steps and
    ## similarly from guess to upper in 100 geometric steps searching
    ## for two points where the function value has opposite signs so
    ## that we can use uniroot to find the root
    left <- (guess / lower) ^ (1 / nstep)
    right <- (upper / guess) ^ (1 / nstep)
    low <- guess
    up <- guess
    L <- NA
    R <- NA
    guess.sign <- sign(f(guess))
    for(i in 1:nstep){
        low <- low / left
        if(sign(f(low)) != guess.sign){
            ## root has been bracketed. Set L & R and exit for loop
            L = low
            R = low * left
            break
        }
        up <- up * right
        if(sign(f(up)) != guess.sign){
            ## root has been bracketed. Set L & R and exit for loop
            L = up / right 
            R = up 
            break 
        }
    }
    if (is.na(L)){
        ## could not bracket root, cannot call uniroot so return NA
        warning("bisection.root failed to locate sign change.")
        return (NA)                         
    }else{
        res <- NA
        tryCatch(
            res <- uniroot(f, c(L, R), tol=toler),
            warning = function(war){}
        )
        if(is.na(res[1])){
            warning("bisection.root failed because uniroot failed.")
            return (NA)
        }
        return (res$root)
    }
}

##' A Newton Raphson root finder: finds x such that f(x) = 0 
##' 
##' The function newton.raphson.root is a general root finder which
##' can find the zero of any function whose derivative is available.
##' In this package, it is called by \code{\link{irr.solve}} and by
##' \code{\link{GenBSImplied}}. It can be used in other situations as
##' well - see the examples below.
##' 
##' @param f The function whose zero is to be found. An R function
##'     object that takes one numeric argument and returns a list of
##'     two components (value and gradient). In an IRR application,
##'     these two components will be the NPV and the DV01/10000. In an
##'     implied volatility application, the components will be the
##'     option price and the vega. See also the examples below
##' @param guess The starting value (guess) from which the solver
##'     starts searching for the IRR
##' @param lower The lower end of the interval within which to search
##'     for the root
##' @param upper The upper end of the interval within which to search
##'     for the root
##' @param max.iter The maximum number of iterations of the
##'     Newton-Raphson procedure
##' @param convergence The relative tolerance threshold used to
##'     determine whether the Newton-Raphson procedure has
##'     converged. The procedure terminates when the last step is less
##'     than \code{convergence} times the current estimate of the
##'     root. Convergence can take place to a non zero local
##'     minimum. This is checked using the \code{toler} criterion
##'     below
##' @param toler The criterion to determine whether a zero has been
##'     found. If the value of the function exceeds \code{toler} in
##'     absolute value, then \code{NA} is returned with a warning
##' @return
##' 
##' The function returns \code{NA} under either of two conditions: (a)
##' the procedure did not converge after \code{max.iter} iterations,
##' or (b) the procedure converged but the function value is not zero
##' within the limits of \code{toler} at this point. The second
##' condition usually implies that the procedure has converged to a
##' non zero local minimum from which there is no downhill gradient.
##'
##' If the iterations converge to a genuine root (within the limits of
##' \code{toler}), then it returns the root that was found.
##'
##' @references The Newton Raphson solver was converted from C++ code
##'     in the \href{http://www.boost.org/}{Boost library}
##' @export
newton.raphson.root <- function(f, guess=0, lower=-Inf, upper=Inf,
                                max.iter=100, toler=1e-6,
                                convergence=1e-8){
    my.warning <- function(reason){
        warning(paste("newton.raphson.root failed:", reason),
                call. = FALSE)
    }
    ## This code is converted from C++ code in Boost Library
    result = guess
    delta = 1
    delta1 = Inf
    delta2 = Inf
    f0 = 0
    iter = 0
    repeat{
        F = f(result)
        f0 = F$value
        f1 = F$gradient
        if(! is.finite(f0) || ! is.finite(f1)){
            my.warning("Function returned non finite value")
            return(NA)
        }
        last.f0 = f0
        delta2 = delta1
        delta1 = delta
        if(0 == f0)
            break
        if(f1 == 0){
            ## Oops zero derivative!!!
            .handle.zero.derivative(f, last.f0, f0, delta, result,
                                    lower, upper)
        }else{
            delta = f0 / f1
        }
        if(abs(delta * 2) > abs(delta2) &&
           is.finite(if (delta > 0) lower else upper)){
            ## last two steps haven't converged, try bisection:
            delta = if (delta > 0){
                (result - lower) / 2
                }else{
                (result - upper) / 2
                }
        }
        guess = result
        result = result - delta  
        if(result <= lower){
            delta = 0.5 * (guess - lower)
            result = guess - delta
            if((result == lower) || (result == upper))
                break
        }else if(result >= upper) {
            delta = 0.5 * (guess - upper)
            result = guess - delta
            if((result == lower) || (result == upper))
                break
        }
        ## update brackets:
        if(delta > 0){
            upper = guess
        }else{
            lower = guess
        }
        iter = iter + 1
        if(iter > max.iter || (abs(delta) < abs(result * convergence)))
            break
    }
    if(abs(f(result)$value) < toler){
        return(result)
    }
    if(iter > max.iter){
        my.warning("Maximum iterations exceeded");
        return(NA)
    }
    if(abs(delta) < abs(result * convergence)){
        my.warning("Converged to a non zero local minimum");
        return(NA)
    }
}

## This function is for internal use by newton.raphson.root and is not
## exported
.handle.zero.derivative <- function(f, last.f0, f0, delta, result,
                                    lower, upper){
    ## This code is converted from C++ code in Boost Library
    if(last.f0 == 0){
        ## this must be the first iteration, pretend that we had a
        ## previous one at either lower or upper:
        if(result == lower){
            guess = upper
        }else{
            guess = lower
        }
        last.f0 = f(guess)$value
        delta = guess - result
    }
    if(sign(last.f0) * sign(f0) < 0){
        ## we've crossed over so move in opposite direction to last
        ## step:
        if(delta < 0){
            delta = (result - lower) / 2
        }else{
            delta = (result - upper) / 2
        }
    }else{
        ## move in same direction as last step:
        if(delta < 0){
            delta = (result - upper) / 2
        }else{
            delta = (result - lower) / 2
        }
   }
}
