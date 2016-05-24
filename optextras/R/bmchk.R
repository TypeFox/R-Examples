bmchk <- function(par, lower = NULL, upper = NULL, 
    bdmsk = NULL, trace = 0, tol = NULL, shift2bound = TRUE) {
    ## Bounds and masks check
    #  ?? check use of par and bvec -- can we simplify?
    # 20101031 -- issue of bounds not working correctly
    #  - -Inf seems to upset bounds
    #  - no upper bounds gives troubles (applies to Rcgmin too!)
    #
    # Input:
    #  par = a vector containing the starting point
    #  lower = vector of lower bounds on parameters
    #  upper = vector of upper bounds on parameters
    # Note: free parameters outside bounds will be adjusted to
    #   bounds.
    # bdmsk = control vector for bounds and masks. Parameters
    #   for which bdmsk are 1
    # are unconstrained or 'free', those with bdmsk 0 are
    #   masked i.e., fixed.
    # For historical reasons, we use the same array as an
    #   indicator that a parameter is at a lower bound (-3)
    #   or upper bound (-1)
    # trace = control of output: 0 for none (default), >0 for
    #   output
    # shift2bound = TRUE if we adjust par values so they are
    #   feasible
    ##
    # Output:
    #    A list with components:
    #     bvec: The parameters adjusted to the nearest bound.
    #     bdmsk: adjusted input masks
    #     bchar: indicator for humans -- "-","L","F","U","+","M"
    #        for out-of-bounds-low, lower bound, free, 
    #            upper bound, out-of-bounds-high, masked (fixed)
    #     lower: adjusted lower bounds
    #     upper: adjusted upper bounds
    #     nolower: TRUE if no lower bounds, FALSE otherwise
    #     noupper: TRUE if no upper bounds, FALSE otherwise
    #     bounds:  TRUE if any bounds, FALSE otherwise
    #     admissible: TRUE if admissible, FALSE if not
    #        No lower bound exceeds an upper bound. That is the 
    #        bounds themselves are sensible. This condition has 
    #        nothing to do with the starting parameters.
    #     maskadded: TRUE if a mask is added, FALSE if not
    #        This implies very close upper and lower bounds for 
    #        the parameters. See the code for the implementation.
    #     parchanged: TRUE if parameters changed, FALSE if not
    #        parchanged = TRUE means that parameters are 
    #        INFEASIBLE, or they would not be changed.
    #     onbound: TRUE if any parameter equal to a bound
    #     feasible: TRUE if feasible, FALSE otherwise
    #
    ########## length of vectors #########
    n <- length(par)
    bvec <- par
    ############# bounds and masks ################
    # set default masks if not defined
    bchar <- rep("F",n) # make sure these are defined
    if (is.null(bdmsk)) {
        bdmsk <- rep(1, n)
    }
    if (trace > 2) {
        cat("bdmsk:")
        print(bdmsk)
    }
    # check if there are bounds
    if (is.null(lower) || !any(is.finite(lower))) 
        nolower <- TRUE
    else nolower <- FALSE
    if (is.null(upper) || !any(is.finite(upper))) 
        noupper <- TRUE
    else noupper <- FALSE
    if (nolower && noupper && all(bdmsk == 1)) 
        bounds <- FALSE
    else bounds <- TRUE
    if (trace > 2) 
        cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, 
            " bounds = ", bounds, "\n")
    if (nolower) 
        lower <- rep(-Inf, n)
    if (noupper) 
        upper <- rep(Inf, n)

## adjust tolerance for masks and parameters ON bounds
    if (is.null(tol)) {
       tol <- .Machine$double.eps * max(abs(par), 1) # use par as bounds Inf
    }
    ######## check bounds and masks #############
    parchanged <- FALSE  ## must be set BEFORE if (bounds) ...; 
                         ## equivalent to feasible<-TRUE
    feasible <- TRUE
    admissible <- TRUE  # similarly must set before we look at bounds
    maskadded <- FALSE  # similarly set here
    onbound <- FALSE
    if (bounds) {
        # Make sure to expand lower and upper
        if (!nolower & (length(lower) < n)) 
            {
                if (length(lower) == 1) {
                  lower <- rep(lower, n)
                }
                else {
                  stop("1<length(lower)<n")
                }
            }  # else lower OK
        if (!noupper & (length(upper) < n)) 
            {
                ## tmp<-readline('Check length upper ')
                if (length(upper) == 1) {
                  upper <- rep(upper, n)
                }
                else {
                  stop("1<length(upper)<n")
                }
            }  # else upper OK
        # At this point, we have full bounds in play
        ######## check admissibility ########
        if (any(lower[which(bdmsk != 0)] > upper[which(bdmsk != 0)])) admissible <- FALSE
        if (trace > 0) cat("admissible = ", admissible, "\n")
        if (any((upper - lower) < tol)) { # essentially masked
            makemask<-which(upper - lower < tol)
            if (trace > 0) {
               cat("Imposing mask as lower ~= upper for following parameters\n")
               print(makemask)
            }
            bdmsk[makemask] <- 0
            bchar[makemask] <- "M"
            # lower[makemask]<- -Inf
            # upper[makemask]<-  Inf
            maskadded <- TRUE
        }
        if (trace > 0) cat("maskadded = ", maskadded, "\n")
        ######## check feasibility ########
        if (admissible) {# This implementation as a loop, but try later to vectorize
            for (i in 1:n) {
                if (bdmsk[i] == 0) {
                  bchar[i] <- "M"
                  # NOTE: we do not change masked parameters, even if out of
                  #   bounds
                  if (!nolower) {
                    if (bvec[i] < lower[i]) {
                      if (trace > 0) {
                        cat("WARNING: ", bvec[i], " = MASKED x[", 
                          i, "] < lower bound = ", lower[i], "\n")
                      }
                      feasible <- FALSE
                    }
                  }
                  if (!noupper) {
                    if (bvec[i] > upper[i]) {
                        if (trace > 0){
                        cat("WARNING: ", bvec[i], " = MASKED x[", 
                        i, "] > upper bound = ", upper[i], "\n")
                        }
                        feasible<-FALSE
                    }
                  }
                }
                else { # not masked, so must be free or active constraint
                  if (!nolower) {
                    if (bvec[i] <= lower[i]) {
                      # Gave trouble 130924 -- <= not < -- bmtest in nlpor
                      # changed 090814 to ensure bdmsk is set; 110105 < not <=
                      if (bvec[i] != lower[i]){
                         bdmsk[i] <- -3.5 # OUT OF BOUNDS LOW
                         bchar[i] <- "-"
                         feasible<-FALSE
                         if (shift2bound) {
                            parchanged <- TRUE
                            if (trace > 0) cat("WARNING: x[", i, "], set ", 
                                 bvec[i], " to lower bound = ", lower[i], "\n")
                            bvec[i] <- lower[i]
                            bdmsk[i] <- -3 
                            bchar[i] <- "L"
                         }
                      } else {
                         bdmsk[i] <- -3  # active lower bound
                         bchar[i] <- "L"
                      }
                    }
                  } # nolower
                  if (!noupper) {
                    if (bvec[i] >= upper[i]) {
                      # changed 090814 to ensure bdmsk is set; 110105 > not >=
                      if (bvec[i] != upper[i]){
                         bdmsk[i] <- -0.5 # OUT OF BOUNDS HIGH
                         bchar[i] <- "+"
                         feasible<-FALSE
                         if (shift2bound) {
                           parchanged <- TRUE
                           if (trace > 0) cat("WARNING: x[", i, "], set ",
                                  bvec[i], " to upper bound = ", upper[i], "\n")
                           bvec[i] <- upper[i]
                           bdmsk[i] <- -1   # active upper bound
                           bchar[i] <- "U"
                         }
                      } else { 
                         bdmsk[i] <- -1  # active upper bound
                         bchar[i] <- "U"
                      }
                    }
                  } # noupper
               }  # end not masked
               ## Should NOT proceed with optimization if feasible and parchanged 
               ## both FALSE (out of bounds)
            }  # end loop for bound/mask check
        } # if admissible
    } # if bounds
    if (trace > 0) 
        cat("parchanged = ", parchanged, "\n")
    if (any(bvec == lower) || any(bvec == upper)){
        onbound <- TRUE
        if (trace > 0) cat("At least one parameter is on a bound\n")
    }
    ############## end bounds check #############
    bcout <- list(bvec, bdmsk, bchar, lower, upper, nolower, noupper, 
        bounds, admissible, maskadded, parchanged, feasible, onbound)
    names(bcout) <- c("bvec", "bdmsk", "bchar", "lower", "upper", "nolower", "noupper", 
       "bounds", "admissible", "maskadded", "parchanged", "feasible", "onbound")
    # Note bdmsk, lower and upper are returned because they are
    #   modified (length, etc.)
    return(bcout)
}  ## end of bmchk.R 
