ugHgenb <- function(par, fnuser=NULL, bdmsk = NULL, lower = NULL, upper = NULL, 
     numgrad=FALSE, control = list()) {
    # Generate the gradient and Hessian for a function specified via the user
    # function wrapper at the parameters par
    #  WITH recognition of bounds and masks.
    #  Uses scaling and maximize
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
    #  bdmsk = index of bounds and masks (see ?? for explanation)
    #  lower = vector of lower bounds on the parameters par
    #  upper = vector of upper bounds on the parameters par
    #  control=list of controls:
    #      asymtol = Tolerance to decide if we should force symmetry
    #                Default is 1.0e-7.
    #      ktrace = integer flag (default 0) to cause allow output if >0
    #      dowarn = TRUE for warnings, FALSE otherwise
    #      stoponerror = logical flag (default=FALSE) to stop if we cannot
    #            compute gradient or Hessian approximations. Otherwise return
    #            with flags gradOK and/or hessOK set FALSE.
    #  ...  = additional arguments to the objective, gradient and Hessian
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
    #  At 2011-10-28 we are NOT using the bounds unless constraints active.
    #
    #################################################################
    if (is.null(fnuser)) stop("MUST have environment fnuser defined.")
    ctrl <- list(asymtol = 1e-07, ktrace = 0, dowarn=TRUE, 
              stoponerror = FALSE)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$ktrace>0) ctrl$dowarn<-TRUE # force true if tracing
    # For bounds constraints, we need to 'project' the gradient and Hessian
    # For masks there is no possibility of movement in parameter.
    bmset <- sort(unique(c(which(par <= lower), which(par >= upper), 
               c(which(bdmsk == 0)))))
    nbm <- length(bmset)  # number of elements nulled by bounds
    #  Note: we assume that we are ON, not over boundary, but use <= and >=.
    # No tolerance is used. ?? Do we want to consider parameters 'close' to
    #   bounds?
    gradOK <- FALSE
    if (ctrl$ktrace > 1) {
        cat("Compute gradient approximation\n")
    }
    if (is.null(fnuser$gr) || numgrad) {
          gt <- try(gn<-grad(ufn, par, fnuser=fnuser), silent = TRUE)
    } else {
          gt <- try(gn<-ugr(par, fnuser=fnuser), silent = TRUE)  
    }
    if (class(gt) != "try-error") gradOK <- TRUE  
    if (!gradOK && ctrl$stoponerror) stop("Gradient computations failure!")
    if (ctrl$ktrace > 0) {
         cat("gradient approx.:\n")
         print(gn)
    }
    if (nbm > 0) {
        gn[bmset] <- 0  # 'project' the gradient
        if (ctrl$ktrace > 0) {
            cat("Adjusted gradient approx.:\n")
            print(gn)
        }
    }
    if (ctrl$ktrace > 0) cat("Compute Hessian approximation\n")
    hessOK <- FALSE
    if (is.null(fnuser$hess)) {
        if (ctrl$ktrace > 1) 
            cat("is.null(hess) is TRUE\n")  # ???
        if (is.null(fnuser$gr) || numgrad) {
            if (ctrl$ktrace > 1) cat("is.null(fnuser$gr) is TRUE use numDeriv hessian()\n")
            Ht<-try(Hn<-hessian(ufn, par, fnuser=fnuser), silent = TRUE) 
            if (class(Ht) == "try-error") {
               if (ctrl$stoponerror) 
                  stop("Unable to compute Hessian using numDeriv::hessian")
               # hessOK still FALSE
            } else hessOK <- TRUE  # Do not need to check for symmetry either.
        } else { # fnuser$gr not null -- use Jacobian
            if (ctrl$ktrace > 1) cat("is.null(fnuser$gr) is FALSE use numDeriv jacobian()\n")
            Ht <- try(Hn <- jacobian(ugr, par, fnuser=fnuser), silent = TRUE)  
            if (class(Ht) == "try-error") {
                if (ctrl$stoponerror) 
                  stop("Unable to compute Hessian using numDeriv::jacobian")
                if (ctrl$ktrace > 0) 
                  cat("Unable to compute Hessian using numDeriv::jacobian \n")
            }
            else hessOK <- TRUE
            if (ctrl$ktrace > 1) {
                cat("Hessian from Jacobian:")
                # print(Hn) # Printed below
            }
            if (!isSymmetric(Hn, tol=sqrt(.Machine$double.eps))) 
                {
                  asym <- sum(abs(t(Hn) - Hn))/sum(abs(Hn))
                  asw <- paste("Hn from jacobian is reported non-symmetric with asymmetry ratio ", 
                    asym, sep = "")
                  if (ctrl$ktrace > 1) 
                    cat(asw, "\n")
                  if (ctrl$dowarn) warning(asw)
                  if (asym > ctrl$asymtol) {
                    if (ctrl$stoponerror) 
                      stop("Hessian too asymmetric")
                  }
                  else hessOK <- TRUE
                  if (ctrl$ktrace > 1) 
                    cat("Force Hessian symmetric\n")
                  else if (ctrl$dowarn) warning("Hessian forced symmetric", call. = FALSE)
                  Hn <- 0.5 * (t(Hn) + Hn)
                }  # end if ! isSymmetric
        }  # numerical hessian at 'solution'
    }
    else { # fnuser$hess exists -- use it
        if (ctrl$ktrace > 1) cat("is.null(fnuser$hess) is FALSE -- trying hess()\n") 
        Hn <- try(uhess(par, fnuser=fnuser), silent = TRUE)  # change 20110222
        if (class(Hn) == "try-error") {
            if (ctrl$stoponerror) 
                stop("Hessian evaluation with function hess() failed")
            if (ctrl$ktrace > 0) 
                cat("Hessian evaluation with function hess() failed \n")
        }
        else {
          # print(Hn)
           if (!isSymmetric(Hn)) {
                asym <- sum(abs(t(Hn) - Hn))/sum(abs(Hn))
                asw <- paste("Hn from hess() is reported non-symmetric with asymmetry ratio ", 
                  asym, sep = "")
                if (ctrl$ktrace > 1) 
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
           else hessOK=TRUE
        }
    }  # end hessian computation
    if (ctrl$ktrace > 1) {
        cat("Raw Hessian with hessOK=",hessOK,"\n")
        print(Hn)
    }
    # apply the bounds and masks
    if (nbm > 0) {
        Hn[bmset, ] <- 0
        Hn[, bmset] <- 0  # and the Hessian
        if (ctrl$ktrace > 1) {
            cat("Adjusted Hessian:\n")
            print(Hn)
        }
    }
    # Scaling SHOULD be taken care of in ufn, ugr, uhess
    ansout <- list(gn, Hn, gradOK, hessOK, nbm)
    names(ansout) <- c("gn", "Hn", "gradOK", "hessOK", "nbm")
    return(ansout)
}  ## end ugHgenb
