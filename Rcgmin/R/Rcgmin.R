Rcgmin <- function(par, fn, gr = NULL, lower = NULL, 
    upper = NULL, bdmsk = NULL, control = list(), ...) {
    ## An R version of the conjugate gradient minimization
    ## using the Dai-Yuan ideas
    #  This version is for unconstrained functions.
    #
    # Input:
    #  par  = a vector containing the starting point
    # fn = objective function (assumed to be sufficeintly
    #   differentiable)
    #  gr = gradient of objective function
    #  lower = vector of lower bounds on parameters
    #  upper = vector of upper bounds on parameters
    # Note: free parameters outside bounds will be adjusted to
    #   bounds.
    # bdmsk = control vector for bounds and masks. Parameters
    #   for which bdmsk are 1
    # are unconstrained or 'free', those with bdmsk 0 are
    #   masked i.e., fixed.
    # For historical reasons, we use the same array as an
    #   indicator that a
    #         parameter is at a lower bound (-3) or upper bound (-1)
    #  control = list of control parameters
    #           maxit = a limit on the number of iterations (default 500)
    #           maximize = TRUE to maximize the function (default FALSE)
    #           trace = 0 (default) for no output,
    #                  >0 for output (bigger => more output)
    # eps=1.0e-7 (default) for use in computing numerical
    #   gradient approximations.
    # dowarn=TRUE by default. Set FALSE to suppress warnings.
    #
    # Output:
    #    A list with components:
    #
    #     par: The best set of parameters found.
    #
    #   value: The value of 'fn' corresponding to 'par'.
    #
    # counts: A two-element integer vector giving the number of
    #   calls to
    # 'fn' and 'gr' respectively. This excludes those calls
    #   needed
    # to compute the Hessian, if requested, and any calls to
    #   'fn'
    # to compute a finite-difference approximation to the
    #   gradient.
    #
    # convergence: An integer code. '0' indicates successful
    #   convergence.
    #          Error codes are
    #          '0' converged
    # '1' indicates that the function evaluation count
    #   'maxfeval'
    #               was reached.
    #          '2' indicates initial point is infeasible
    #
    # message: A character string giving any additional
    #   information returned
    #          by the optimizer, or 'NULL'.
    #
    # bdmsk: Returned index describing the status of bounds and
    #   masks at the
    # proposed solution. Parameters for which bdmsk are 1 are
    #   unconstrained
    # or 'free', those with bdmsk 0 are masked i.e., fixed. For
    #   historical
    # reasons, we indicate a parameter is at a lower bound
    #   using -3
    #          or upper bound using -1.
    #
    #
    #  Author:  John C Nash
    #  Date:  April 2, 2009; revised July 28, 2009
    #################################################################
    # control defaults -- idea from spg
    if (is.null(control$trace)) control$trace=0
    # check if there are bounds
    if (is.null(lower) || !any(is.finite(lower))) 
        nolower = TRUE
    else nolower = FALSE
    if (is.null(upper) || !any(is.finite(upper))) 
        noupper = TRUE
    else noupper = FALSE
    if (nolower && noupper && all(bdmsk == 1)) 
        bounds = FALSE
    else bounds = TRUE
    if (control$trace > 1) 
        cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, 
            " bounds = ", bounds, "\n")
    grNULL <- is.null(gr)
    if (is.null(control$dowarn)) control$dowarn<-TRUE
    #############################################
    if (bounds) {
       ans<-Rcgminb(par, fn, gr, lower = lower, 
          upper = upper, bdmsk = bdmsk, control = control, ...)
    } else {
       ans<-Rcgminu(par, fn, gr, control = control, ...)
    }
    return(ans) # ?? is this needed
}  ## end of Rcgmin
