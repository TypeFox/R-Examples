gHgen <- function(par, fn, gr = NULL, hess = NULL, control = list(ktrace=0), 
    ...) {
    # Generate the gradient and Hessian for a given function at the parameters
    #   par.
    #
    # Input:
    #  par = a vector for function parameters
    #  fn = a user function (assumed to be sufficeintly differentiable)
    # gr = name of a function to compute the (analytic) gradient of the user
    #   function
    # hess = name of a function to compute the (analytic) hessian of the user
    #   function
    #         This will rarely be available, but is included for completeness.
    #  control=list of controls:
    #      asymtol = Tolerance to decide if we should force symmetry
    #                Default is 1.0e-7.
    #      ktrace = logical flag (default FALSE) to monitor progress
    #      stoponerror = logical flag (default=FALSE) to stop if we cannot
    #            compute gradient or Hessian approximations. Otherwise return
    #            with flags gradOK and/or hessOK set FALSE.
    # ...  = additional arguments to the objective, gradient and Hessian
    #   functions
    #
    # Output:
    # A list containing
    #    g = a vector containing the gradient
    #    H = a matrix containing the Hessian
    #    gradOK = logical indicator that gradient is OK. TRUE/FALSE
    #    hessOK = logical indicator that Hessian is OK. TRUE/FALSE
    #    nbm = 0: Number of bounds and masks active. Not used in gHgen.
    #
    #  Author: John Nash
    #  Date:  January 14, 2011; updated June 25, 2011
    #
    #################################################################
    ctrl <- list(asymtol = 1e-07, ktrace = FALSE, stoponerror = FALSE)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    gradOK <- FALSE
    if (ctrl$ktrace) 
        cat("Compute gradient approximation\n")
    if (is.null(gr)) {
        gn <- try(grad(fn, par, ...), silent = TRUE)  # change 20100711
    }
    else {
        gn <- try(gr(par, ...), silent = TRUE)  # Gradient at solution # change 20100711
    }
    if (class(gn) != "try-error") 
        gradOK <- TRUE  # 100215 had == rather than != here
    if (!gradOK && ctrl$stoponerror) 
        stop("Gradient computations failure!")
    if (ctrl$ktrace) 
        print(gn)
    if (ctrl$ktrace) 
        cat("Compute Hessian approximation\n")
    hessOK <- FALSE
    if (is.null(hess)) {
        if (is.null(gr)) {
            if (ctrl$ktrace>0) cat("is.null(gr) is TRUE use numDeriv hessian()\n") # ???
            Hn <- try(hessian(fn, par, ...), silent = TRUE)  # change 20100711
            if (class(Hn) == "try-error") {
                if (ctrl$stoponerror) 
                  stop("Unable to compute Hessian using numDeriv::hessian")
                # hessOK still FALSE
            }
            else hessOK <- TRUE  # Do not need to check for symmetry either.
        }
        else {
            if (ctrl$ktrace>0) cat("is.null(gr) is FALSE use numDeriv jacobian()\n") # ???
            Hn <- try(jacobian(gr, par, ...), silent = TRUE)  # change 20100711
            if (class(Hn) == "try-error") {
                if (ctrl$stoponerror) 
                  stop("Unable to compute Hessian using numderiv::jacobian")
            }
            else hessOK <- TRUE
            if (ctrl$ktrace>0) {
               cat("Hessian from jacobian:")
               print(Hn)
            }
            if (!isSymmetric(Hn)) 
                {
                  asym <- sum(abs(t(Hn) - Hn))/sum(abs(Hn))
                  asw <- paste("Hn from jacobian is reported non-symmetric with asymmetry ratio ", 
                    asym, sep = "")
                  if (ctrl$ktrace) 
                    cat(asw, "\n")
                  warning(asw)
                  cat("asym, ctrl$asymtol: ", asym, ctrl$asymtol, "\n")
                  if (asym > ctrl$asymtol) {
                    if (ctrl$stoponerror) 
                      stop("Hessian too asymmetric")
                  }
                  else hessOK <- TRUE
                  if (ctrl$ktrace) 
                    cat("Force Hessian symmetric\n")
                  else warning("Hessian forced symmetric")
                  Hn <- 0.5 * (t(Hn) + Hn)
                }  # end if ! isSymmetric
        }  # numerical hessian at 'solution'
    }
    else {
        if (ctrl$ktrace>0) cat("is.null(hess) is FALSE -- trying hess()\n") # ???
        Hn <- try(hess(par, ...), silent = TRUE)  # change 20110222
        if (class(Hn) == "try-error") {
            if (ctrl$stoponerror) 
                stop("Hessian evaluation with function hess() failed")
        }
        else hessOK <- TRUE
        if (!isSymmetric(Hn, tol=0.01*sqrt(.Machine$double.eps))) 
            {
                asym <- sum(abs(t(Hn) - Hn))/sum(abs(Hn))
                asw <- paste("Hn from hess() is reported non-symmetric with asymmetry ratio ", 
                  asym, sep = "")
                if (ctrl$ktrace) 
                  cat(asw, "\n")
                warning(asw)
                if (asym > ctrl$asymtol) {
                  if (ctrl$stoponerror) 
                    stop("Hessian too asymmetric")
                }
                else hessOK <- TRUE
                if (ctrl$ktrace) 
                  cat("Force Hessian symmetric\n")
                else warning("Hessian forced symmetric")
                Hn <- 0.5 * (t(Hn) + Hn)
            }  # end if ! isSymmetric
    }  # end hessian computation
    if (ctrl$ktrace) 
        print(Hn)
    ansout <- list(gn, Hn, gradOK, hessOK, 0)
    names(ansout) <- c("gn", "Hn", "gradOK", "hessOK", "nbm")
    return(ansout)
}  ## end gHgen
