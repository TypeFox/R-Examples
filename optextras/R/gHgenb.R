gHgenb <- function(par, fn, gr = NULL, hess = NULL, bdmsk = NULL, 
    lower = NULL, upper = NULL, control = list(ktrace=0), ...) {
    # Generate the gradient and Hessian for a given function at the parameters
    #   par
    #  WITH recognition of bounds and masks.
    #################################################################
    #  NOTE: This routine generates the gradient and Hessian ignoring
    #    bounds and masks and THEN applies a projection. This can cause
    #    difficulties if the small steps taken to compute derivatives
    #    give inadmissible parameter sets to functions.
    #    We also apply symmetry check BEFORE the bounds/masks projection.
    #    JN 2011-6-25
    #    ?? How should that be dealt with generally?? 110625
    #    ?? specifically can do 1 parameter at a time
    #################################################################
    #
    # Input:
    #  par = a vector for function parameters
    #  fn = a user function (assumed to be sufficeintly differentiable)
    # gr = name of a function to compute the (analytic) gradient of the user
    #   function
    # hess = name of a function to compute the (analytic) hessian of the user
    #   function
    #         This will rarely be available, but is included for completeness.
    #  bdmsk = index of bounds and masks (see ?? for explanation)
    #  lower = vector of lower bounds on the parameters par
    #  upper = vector of upper bounds on the parameters par
    #  control=list of controls:
    #      asymtol = Tolerance to decide if we should force symmetry
    #                Default is 1.0e-7.
    #      ktrace = integer flag (default 0) to cause allow output if >0
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
    #
    #  Author: John Nash
    #  Date:  January 14, 2011; updated June 25, 2011
    #
    #################################################################
    #  IMPORTANT: There is a serious question of how to deal with Hessian
    #    and gradient at active bounds, since these have meaning on one
    #    side of the constraint only. JN 2011-6-25
    #
    #################################################################
    ctrl <- list(asymtol = 1e-07, ktrace = 0, dowarn=TRUE, stoponerror = FALSE)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$ktrace>0) ctrl$dowarn<-TRUE # force TRUE for trace
    # For bounds constraints, we need to 'project' the gradient and Hessian
    # For masks there is no possibility of movement in parameter.
    bmset <- sort(unique(c(which(par <= lower), which(par >= upper), c(which(bdmsk == 
        0)))))
    nbm <- length(bmset)  # number of elements nulled by bounds
    #  Note: we assume that we are ON, not over boundary, but use <= and >=.
    # No tolerance is used. ?? Do we want to consider parameters 'close' to
    #   bounds?
    gradOK <- FALSE
    if (ctrl$ktrace > 0) 
        cat("Compute gradient approximation\n")
    if (is.null(gr)) {
        gn <- try(grad(fn, par, ...), silent = TRUE)  # change 20100711
    }
    else {
        gn <- try(gr(par, ...), silent = TRUE)  # Gradient at solution # change 20100711
    }
    if (class(gn) != "try-error") 
        gradOK <- TRUE  
    if (!gradOK && ctrl$stoponerror) 
        stop("Gradient computations failure!")
    if (ctrl$ktrace > 0)  
        print(gn)
    if (nbm > 0) {
        gn[bmset] <- 0  # 'project' the gradient
        if (ctrl$ktrace > 0) {
            cat("Adjusted gradient:\n")
            print(gn)
        }
    }
    if (ctrl$ktrace > 0) 
        cat("Compute Hessian approximation\n")
    hessOK <- FALSE
    if (is.null(hess)) {
        if (ctrl$ktrace > 0) 
            cat("is.null(hess) is TRUE\n")  # ???
        if (is.null(gr)) {
            if (ctrl$ktrace > 0) 
                cat("is.null(gr) is TRUE use numDeriv hessian()\n")  # ???
            Hn <- try(hessian(fn, par, ...), silent = TRUE)  # change 20100711
            if (class(Hn) == "try-error") {
                if (ctrl$stoponerror) 
                  stop("Unable to compute Hessian using numDeriv::hessian")
                # hessOK still FALSE
            }
            else hessOK <- TRUE  # Do not need to check for symmetry either.
        }
        else {
            if (ctrl$ktrace > 0) 
                cat("is.null(gr) is FALSE use numDeriv jacobian()\n")  # ???
            tHn <- try(Hn <- jacobian(gr, par, ...), silent = TRUE)  # change 20100711
            if (class(tHn) == "try-error") {
                if (ctrl$stoponerror) 
                  stop("Unable to compute Hessian using numderiv::jacobian")
                if (ctrl$ktrace > 0) 
                  cat("Unable to compute Hessian using numderiv::jacobian \n")
            }
            else hessOK <- TRUE
            if (ctrl$ktrace > 0) {
                cat("Hessian from Jacobian:")
                # print(Hn) # Printed below
            }
            if (!isSymmetric(Hn, tol=0.1*sqrt(.Machine$double.eps))) 
                {
                  asym <- sum(abs(t(Hn) - Hn))/sum(abs(Hn))
                  asw <- paste("Hn from jacobian is reported non-symmetric with asymmetry ratio ", 
                    asym, sep = "")
                  if (ctrl$ktrace > 0) 
                    cat(asw, "\n")
                  if (ctrl$dowarn) warning(asw)
                  if (asym > ctrl$asymtol) {
                    if (ctrl$stoponerror) 
                      stop("Hessian too asymmetric")
                  }
                  else hessOK <- TRUE
                  if (ctrl$ktrace > 0) 
                    cat("Force Hessian symmetric\n")
                  else if (ctrl$dowarn) warning("Hessian forced symmetric", call. = FALSE)
                  Hn <- 0.5 * (t(Hn) + Hn)
                }  # end if ! isSymmetric
        }  # numerical hessian at 'solution'
    }
    else {
        if (ctrl$ktrace > 0) 
            cat("is.null(hess) is FALSE -- trying hess()\n")  # ???
        Hn <- try(hess(par, ...), silent = TRUE)  # change 20110222
        if (class(Hn) == "try-error") {
            if (ctrl$stoponerror) 
                stop("Hessian evaluation with function hess() failed")
            if (ctrl$ktrace > 0) 
                cat("Hessian evaluation with function hess() failed \n")
        }
        else hessOK <- TRUE
        print(Hn)
        if (!isSymmetric(Hn)) 
            {
                asym <- sum(abs(t(Hn) - Hn))/sum(abs(Hn))
                asw <- paste("Hn from hess() is reported non-symmetric with asymmetry ratio ", 
                  asym, sep = "")
                if (ctrl$ktrace > 0) 
                  cat(asw, "\n")
                else if (ctrl$dowarn) warning(asw)
                if (asym > ctrl$asymtol) {
                  if (ctrl$stoponerror) 
                    stop("Hessian too asymmetric")
                }
                else hessOK <- TRUE
                if (ctrl$ktrace > 0) 
                  cat("Force Hessian symmetric\n")
                else if (ctrl$dowarn) warning("Hessian forced symmetric", call. = FALSE)
                Hn <- 0.5 * (t(Hn) + Hn)
            }  # end if ! isSymmetric
    }  # end hessian computation
    if (ctrl$ktrace > 0) 
        print(Hn)
    # apply the bounds and masks
    if (nbm > 0) {
        Hn[bmset, ] <- 0
        Hn[, bmset] <- 0  # and the Hessian
        if (ctrl$ktrace > 0) {
            cat("Adjusted Hessian:\n")
            print(Hn)
        }
    }
    ansout <- list(gn, Hn, gradOK, hessOK, nbm)
    names(ansout) <- c("gn", "Hn", "gradOK", "hessOK", "nbm")
    return(ansout)
}  ## end gHgenb
