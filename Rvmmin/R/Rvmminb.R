Rvmminb <- function(par, fn, gr=NULL, lower = NULL, 
    upper = NULL, bdmsk, control = list(), ...) {
    #1
    ## Bounds constrained version -- 090408, 090814
    # 20101031 -- issue of bounds not working correctly
    #  - -Inf seems to upset bounds
    #  - no upper bounds gives troubles (applies to Rcgmin too!)
    #
    ## An R version of the Nash version of Fletcher's Variable
    #   Metric minimization
    # This uses a simple backtracking line search.
    # Input:
    #  par = a vector containing the starting point
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
    # maxit = a limit on the gradient evaluations (default
    #   1500)
    # ?? do we want to make it vary with n??
    # maxfeval = a limit on the function evaluations (default
    #   10000)
    #           maximize = TRUE to maximize the function (default FALSE)
    #           trace = 0 (default) for no output,
    #                  >0 for output (bigger => more output)
    #
    # dowarn=TRUE by default. Set FALSE to suppress warnings.
    #
    ##
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
    #
    # '1' indicates that the iteration limit 'maxit' or the
    #   function
    #               evaluation limit mafeval have been reached.
    # '20' indicates inadmissible input parameters for function
    #   (bounds may be OK)
    # '21' indicates inadmissible intermediate parameters for
    #   function
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
    #  Author:  John C Nash
    #  Date:  May 14, 2012
    #################################################################
    # control defaults
    # NOT yet in control set ??    #  ?? put keepinputpar into controls??
    ctrl <- list(maxit = 500, maxfeval = 3000, maximize = FALSE, 
        trace = 0, eps = 1e-07, dowarn = TRUE, keepinputpar = FALSE, acctol = 0.0001)
    # keepinputpar TRUE => Do not let bmchk change parameters to nearest bound
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control  #
    maxit <- ctrl$maxit  #
    maxfeval <- ctrl$maxfeval  #
    maximize <- ctrl$maximize  # TRUE to maximize the function
    trace <- ctrl$trace  #
    eps <- ctrl$eps  #
    acctol <- ctrl$acctol
    dowarn <- ctrl$dowarn  #
    fargs <- list(...)  # the ... arguments that are extra function / gradient data
    #################################################################
    ## Set working parameters (See CNM Alg 22)
    if (trace > 0) 
        cat("Rvmminb -- J C Nash 2012 - Box-constrained R version of Alg 21\n")
    bvec <- par  # copy the parameter vector
    n <- length(bvec)  # number of elements in par vector
    if (trace > 0) {
        cat("Problem of size n=", n, "  Dot arguments:\n")
        print(fargs)
    }
    ifn <- 1  # count function evaluations
    stepredn <- 0.2  # Step reduction in line search
    reltest <- 100  # relative equality test
    ceps <- .Machine$double.eps * reltest
    dblmax <- .Machine$double.xmax  # used to flag bad function
    #############################################
    # gr MUST be provided
    if (is.null(gr)) {  # if gr function is not provided STOP (Rvmmin has definition)
       stop("A gradient calculation (analytic or numerical) MUST be provided for Rvmminb") 
    }
    if ( is.character(gr) ) {
       # Convert string to function call, assuming it is a numerical gradient function
       mygr<-function(par=par, userfn=fn, ...){
           do.call(gr, list(par, userfn, ...))
       }
    } else { mygr<-gr }
    ############# end test gr ####################
    f<-try(fn(bvec, ...), silent=TRUE) # Compute the function.
    if ((class(f) == "try-error") | is.na(f) | is.null(f) | is.infinite(f)) {
        msg <- "Initial point gives inadmissible function value"
        conv <- 20
        if (trace > 0) 
            cat(msg, "\n") # change NA to dblmax 110524
        ans <- list(bvec, dblmax, c(ifn, 0), conv, msg, bdmsk)  #
        names(ans) <- c("par", "value", "counts", "convergence", 
            "message", "bdmsk")
        return(ans)
    }
    if (maximize) f <- -f
    if (trace > 0) cat("Initial fn=", f, "\n")
    if (trace > 2) print(bvec)
    keepgoing <- TRUE  # to ensure loop continues until we are finished
    ig <- 1  # count gradient evaluations
    ilast <- ig  # last time we used gradient as search direction
    fmin <- f  # needed for numerical gradients
    g <- mygr(bvec, ...)  # Do we need to use try() ?? Possible not
    if (maximize) g <- -g
    if (trace > 2) {
        cat("g:")
        print(g)
    }
    oldstep <- 1
    conv <- -1
    while (keepgoing) { ## main loop -- must remember to break out of it!
        if (ilast == ig) { # reset the approx. inverse hessian B to unit matrix
            B <- diag(1, n, n)  # create unit matrix of order n
            if (trace > 2) cat("Reset Inv. Hessian approx at ilast = ", ilast, "\n")
        }
        fmin <- f
        if (trace > 0) cat(" ", ifn, " ", ig, " ", fmin, "\n")
        par <- bvec  # save parameters
        if (!all(is.numeric(g))) {
            g <- rep(0, n)  # 110619
            cat("zeroing gradient because of failure\n")
        }
        c <- g  # save gradient
        ## Bounds and masks adjustment of gradient ##
        ## current version with looping -- later try to vectorize
        ##         if (bounds) 
        if (trace > 3) {
             cat("bdmsk:")
             print(bdmsk)
        }
        for (i in 1:n) {
            if ((bdmsk[i] == 0)) {
               g[i] <- 0
            }
            else {
               if (bdmsk[i] == 1) {
                  if (trace > 2) 
                     cat("Parameter ", i, " is free\n")
               }
               else {
                  if ((bdmsk[i] + 2) * g[i] < 0) {
                     g[i] <- 0  # active mask or constraint
                  }
                  else {
                     bdmsk[i] <- 1  # freeing parameter i
                     if (trace > 1) 
                       cat("freeing parameter ", i, "\n")
                  }
               }
            }
        }  # end masking loop on i
               if (trace > 3) {
                  cat("bdmsk adj:")
                  print(bdmsk)
                  cat("proj-g:")
                  print(g)
               }
               ## end bounds and masks adjustment of gradient
        ###    }  # if bounds
        t <- as.vector(-B %*% g)  # compute search direction
        if (!all(is.numeric(t))) 
            t <- rep(0, n)  # 110619
        if (trace > 2) {
            cat("t:")
            print(t)
        }
        t[which(bdmsk <= 0)] <- 0  # apply masks and box constraints
        if (trace > 2) {
            cat("adj-t:")
            print(t)
        }
        gradproj <- sum(t * g)  # gradient projection
        if (trace > 1) 
            cat("Gradproj =", gradproj, "\n")
        accpoint <- FALSE  # Need this BEFORE gradproj test
        if (is.nan(gradproj)) {
            warning("gradproj Nan")
            gradproj <- 0  # force null
        }
        if (gradproj < 0) 
            {
                # Must be going downhill
                ########################################################
                ####      Backtrack only Line search                ####
                changed <- TRUE  # Need to set so loop will start
                steplength <- oldstep
                while ((f >= fmin) && changed && (!accpoint)) {
                   # We seek a lower point, but must change parameters too
                   ###if (bounds) {
                   # Box constraint -- adjust step length for free parameters
                   for (i in 1:n) {
                      # loop on parameters -- vectorize??
                      if ((bdmsk[i] == 1) && (t[i] != 0)) {
                         # only concerned with free parameters and non-zero search
                         #   dimension
                         if (t[i] < 0) {
                            # going down. Look at lower bound
                            trystep <- (lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
                         }
                         else {
                            # going up, check upper bound
                            trystep <- (upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
                         }
                         if (trace > 2) 
                            cat("steplength, trystep:", steplength, trystep, "\n")
                         steplength <- min(steplength, trystep)  # reduce as necessary
                       }  # end steplength reduction
                   }  # end loop on i to reduce step length
                   # end box constraint adjustment of step length
                   if (trace > 1) cat("reset steplength=", steplength, "\n")
                   ###  }  # end if bounds
                   # end box constraint adjustment of step length
                   bvec <- par + steplength * t
                   if (trace > 2) {
                     cat("new bvec:")
                     print(bvec)
                   }
                   changed <- (!identical((bvec + reltest), (par + 
                     reltest)))
                   if (changed) {
                     # compute new step, if possible
                     f <- fn(bvec, ...)  # Because we need the value for linesearch, don't use try()
                    # instead preferring to fail out, which will hopefully be unlikely.
                    if (maximize) f <- -f
                    ifn <- ifn + 1
                    if (ifn > maxfeval) {
                      msg <- "Too many function evaluations"
                      if (dowarn) 
                        warning(msg)
                      conv <- 1
                      changed <- FALSE
                      keepgoing <- FALSE
                      break
                    }
                    if (is.na(f) | is.null(f) | is.infinite(f)) {
                      if (trace > 2) {
                        cat("Function is not calculable at intermediate bvec:")
                        print(bvec)
                      }
                      #               msg='Function is not calculable at an intermediate point'
                      #               #  stop('f is NA')
                      #               conv<-21
                      #               break
                      f <- dblmax  # try big function to escape
                    }
                    if (f < fmin) {
                      # We have a lower point. Is it 'low enough' i.e.,
                      #   acceptable
                      accpoint <- (f <= fmin + gradproj * steplength * 
                        acctol)
                    }
                    else {
                      steplength <- steplength * stepredn
                      if (trace > 0) 
                        cat("*")
                    }
                  }
                  else {
                    # NOT changed in step reduction
                    if (trace > 1) 
                      cat("Unchanged in step redn \n")
                  }
                }  # end while ((f >= fmin) && changed )
            }  # end if gradproj<0
        if (accpoint) {
            # matrix update if acceptable point.
            ### if (bounds) {
            for (i in 1:n) { ## Reactivate constraints??
                if (bdmsk[i] == 1) {
                    # only interested in free parameters
                    # make sure < not <= below to avoid Inf comparisons
                    if ((bvec[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) {
                       # are we near or lower than lower bd
                       if (trace > 2) cat("(re)activate lower bd ", i, " at ", lower[i], "\n")
                       bdmsk[i] <- -3
                    }  # end lower bd reactivate
                    if ((upper[i] - bvec[i]) < ceps * (abs(upper[i]) + 1)) {
                        # are we near or above upper bd
                        if (trace > 2) cat("(re)activate upper bd ", i," at ", upper[i], "\n")
                        bdmsk[i] <- -1
                    }  # end lower bd reactivate
                }  # end test on free params
             }  # end reactivate constraints
             ###   }  # if bounds
            test <- try(g <- mygr(bvec, ...), silent = TRUE)  # ?? use try()
            if (class(test) == "try-error") 
                stop("Bad gradient!!")
            if (any(is.nan(g))) 
                stop("NaN in gradient")
            ig <- ig + 1
            if (maximize) g <- -g
            if (ig > maxit) {
                keepgoing = FALSE
                msg = "Too many gradient evaluations"
                if (dowarn) 
                  warning(msg)
                conv <- 1
                break
            }
            g[which(bdmsk == 0)] <- 0  # active mask or constraint
            t <- as.vector(steplength * t)
            c <- as.vector(g - c)
            D1 <- sum(t * c)
            if (D1 > 0) {
                y <- as.vector(crossprod(B, c))
                D2 <- as.double(1+crossprod(c,y)/D1)  # as.double because D2 is a 1 by 1 matrix otherwise
                # May be able to be more efficient below -- need to use
                #   outer function
                B <- B - (outer(t, y) + outer(y, t) - D2 * outer(t, t))/D1
            }
            else {
                if (trace > 0) 
                  cat("UPDATE NOT POSSIBLE\n")
                ilast <- ig  # note gradient evaluation when update failed
            }  # D1 > 0 test
        }
        else {
            # no acceptable point
            if (trace > 0) 
                cat("No acceptable point\n")
            if (ig == ilast) {
                # we reset to gradient and did new linesearch
                keepgoing <- FALSE  # no progress possible
                if (conv < 0) {
                  conv <- 0
                  msg <- "Converged"
                }
                if (trace > 0) 
                  cat(msg, "\n")
            }
            else {
                ilast <- ig  # reset to gradient search
                if (trace > 0) 
                  cat("Reset to gradient search\n")
            }  # end else ig != ilast
        }  # end else no accpoint
    }  # end main loop  (while keepgoing)
    if (trace > 0) 
        cat("Seem to be done VM\n")
    if (maximize) 
        fmin <- (-1) * fmin
    msg <- "Rvmminb appears to have converged"
    ans <- list(par, fmin, c(ifn, ig), convergence=conv, msg, bdmsk)
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message", "bdmsk")
    #return(ans)
    ans
}  ## end of Rvmminb
