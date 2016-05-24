Rvmminu <- function(par, fn, gr=NULL, control = list(), ...) {
    ## Unconstrained version -- 120521
    ## An R version of the Nash version of Fletcher's Variable
    #   Metric minimization
    # See comments in vm
    #  Author:  John C Nash
    #  Date:  May 21, 2012
    #################################################################
    # control defaults
    # NOT yet in control set ??    #  ?? put keepinputpar into controls??
    ctrl <- list(maxit = 500, maxfeval = 3000, maximize = FALSE, 
        trace = 0, eps = 1e-07, dowarn = TRUE, acctol = 0.0001, checkgrad=TRUE)
    # checkgrad not needed, but here to avoid error
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control  #
    maxit <- ctrl$maxit  #
    maxfeval <- ctrl$maxfeval  #
    maximize <- ctrl$maximize  # TRUE to maximize the function
    trace <- ctrl$trace  #
    eps <- ctrl$eps  #
    acctol <- ctrl$acctol # 130125
    dowarn <- ctrl$dowarn  #
    fargs <- list(...)  # the ... arguments that are extra function / gradient data
    #################################################################
    ## Set working parameters (See CNM Alg 22)
    if (trace > 0) 
        cat("Rvmminu -- J C Nash 2009, 2011, 2012 - an R implementation of Alg 21\n")
    bvec <- par  # copy the parameter vector
    n <- length(bvec)  # number of elements in par vector
    if (trace > 0) {
        cat("Problem of size n=", n, "  Dot arguments:\n")
        print(fargs)
    }
    ifn <- 1  # count function evaluations
    stepredn <- 0.2  # Step reduction in line search
    acctol <- 1e-04  # acceptable point tolerance
    reltest <- 100  # relative equality test
    ceps <- .Machine$double.eps * reltest
    dblmax <- .Machine$double.xmax  # used to flag bad function
    #############################################
    # gr MUST be provided
    if (is.null(gr)) {  # if gr function is not provided STOP (vm has definition)
       stop("A gradient calculation (analytic or numerical) MUST be provided for vmb") 
    }
    if ( is.character(gr) ) {
       # Convert string to function call, assuming it is a numerical gradient function
       mygr<-function(par=par, userfn=fn, ...){
           do.call(gr, list(par, userfn, ...))
       }
    } else { mygr<-gr }
#    cat(deparse(substitute(mygr)),"\n")
#    tmp<-readline("mygr:")
#    print(mygr)
    ############# end test gr ####################
    f<-try(fn(bvec, ...), silent=TRUE) # Compute the function.
    if ((class(f) == "try-error") | is.na(f) | is.null(f) | is.infinite(f)) {
        msg <- "Initial point gives inadmissible function value"
        conv <- 20
        if (trace > 0) 
            cat(msg, "\n") # change NA to dblmax 110524
        ans <- list(bvec, dblmax, c(ifn, 0), 0, conv, msg)  #
        names(ans) <- c("par", "value", "counts", "convergence", 
            "message")
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
        t <- as.vector(-B %*% g)  # compute search direction
        if (!all(is.numeric(t))) 
            t <- rep(0, n)  # 110619
        if (trace > 2) {
            cat("t:")
            print(t)
        }
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
                accpoint <- (f <= fmin + gradproj * steplength * acctol)
                while (changed && (!accpoint)) {
                  # We seek a lower point, but must change parameters too
                  # end box constraint adjustment of step length
                  bvec <- par + steplength * t
                  if (trace > 2) {
                    cat("new bvec:")
                    print(bvec)
                  }
                  changed <- (!identical((bvec + reltest), (par + reltest)))
                  if (changed) {
                    # compute new step, if possible
                    f <- fn(bvec, ...)  # Because we need the value for linesearch, don't use try()
                    # instead preferring to fail out, which will hopefully be unlikely.
                    if (maximize) f <- -f
                    if (trace > 2) {
                       cat("New f=",f,"\n")
                    }
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
#                    tmp<-readline("Test the function")
                    accpoint <- (f <= fmin + gradproj * steplength * acctol)
                    if (! accpoint) {
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
    msg <- "Rvmminu appears to have converged"
    ans <- list(par, fmin, c(ifn, ig), convergence=conv, msg)
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message")
    #return(ans)
    ans
}  ## end of vmu
