Rvmmin <- function(par, fn, gr = NULL, lower = NULL, 
    upper = NULL, bdmsk = NULL, control = list(), ...) {
    ## An R version of the Nash version of Fletcher's Variable
    #   Metric minimization
    # This uses a simple backtracking line search.
    # This is a driver for both constrained and unconstrained routines.
    #
    # Input:
    # par  = a vector containing the starting point
    # fn = objective function (assumed to be sufficeintly
    #   differentiable)
    # gr = gradient of objective function
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
    #           dowarn=TRUE by default. Set FALSE to suppress warnings.
    #           checkgrad = TRUE by default. Check analytic gradient 
    #                  against numDeriv results.
    #           checkbounds = TRUE by default. Check parameters and bounds
    #                  for addmissible bounds and feasible start.
    #
    # Output:
    #    A list with components:
    #
    #   par: The best set of parameters found.
    #
    #   value: The value of 'fn' corresponding to 'par'.
    #
    #   counts: A two-element integer vector giving the number of
    #     calls to
    #     'fn' and 'gr' respectively. This excludes those calls
    #     needed to compute the Hessian, if requested, and any 
    #     calls to 'fn' to compute a finite-difference approximation 
    #     to the gradient.
    #
    #   convergence: An integer code. '0' indicates successful
    #     convergence.
    #          Error codes are
    #          '0' converged
    #          '1' indicates that the function evaluation count
    #              'maxfeval' was reached.
    #          '2' indicates initial point is infeasible
    #
    #   message: A character string giving any additional
    #      information returned by the optimizer, or 'NULL'.
    #
    #   bdmsk: Returned index describing the status of bounds and
    #      masks at the proposed solution. Parameters for which 
    #      bdmsk are 1 are unconstrained or 'free', those with 
    #      bdmsk 0 are masked i.e., fixed. For historical
    #      reasons, we indicate a parameter is at a lower bound
    #      using -3 or upper bound using -1.
    #
    #
    #  Author:  John C Nash
    #  Date:  April 2, 2009; revised July 28, 2009, May 21, 2012
    #################################################################
    #
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
    if (is.null(gr)) {
       gr <- "grfwd" # use forward gradient approximation if no gradient code provided
       if (control$trace > 0) cat("WARNING: forward gradient approximation being used\n")
    } else {
       if (is.character(gr)) { # assume numerical gradient
           if (control$trace > 0) cat("WARNING: using gradient approximation '",gr,"'\n")
       } else { # analytic gradient, so check if requested
           if (is.null(control$checkgrad)) control$checkgrad <- TRUE
           if (control$checkgrad) { # check gradient
              testgrad<-grchk(par, fn, gr, trace=control$trace, ...)
              if (! testgrad) warning("Gradient code for Rvmmin may be faulty - check it!")
           }
       } # end else
    }
    control$checkgrad<-NULL # to avoid problems in subsidiary routines
    if (is.null(control$dowarn)) control$dowarn<-TRUE
    #############################################
    if (bounds) { 
       if (is.null(control$checkbounds)) { control$checkbounds <- TRUE }
       ### Check bounds feasible
       if (control$checkbounds) {
          btest <- bmchk(par, lower = lower, upper = upper, bdmsk = bdmsk, 
             trace = control$trace)
          if (!btest$admissible) 
             stop("Inadmissible bounds: one or more lower>upper")
          if (btest$parchanged) {
             if (is.null(control$keepinputpar) || ! control$keepinputpar) { 
                 warning("Parameter out of bounds has been moved to nearest bound")
                 control$keepinputpar<-NULL # avoid problems in subsidiary routines
             } else stop("Parameter out of bounds")
          }
       }
       nolower <- btest$nolower
       noupper <- btest$noupper
       bounds <- btest$bounds
       bdmsk <- btest$bdmsk  # change bdmsk to values set in bmchk
       if (control$trace > 3) {
          cat("Adjusted bdmsk vector:")
          print(bdmsk)
       }
       lower <- btest$lower
       upper <- btest$upper
       control$checkbounds<-NULL # to avoid problems in subsidiary routines
       ############## end bounds check #############
       ans<-Rvmminb(par, fn, gr, lower = lower, 
          upper = upper, bdmsk = bdmsk, control = control, ...)
    } else {
       ans<-Rvmminu(par, fn, gr, control = control, ...)
    }
#    return(ans) # commented 130108
}  ## end of Rvmmin
