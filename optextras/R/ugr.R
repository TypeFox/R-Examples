############### ugr.R ####################
ugr <- function(par, fnuser) {
    # gradient wrapper. par are INTERNAL scaled parameters == upar/parscale
    npar <- length(par)
    parps<-par*fnuser$PARSCALE # user parameters
#    cat("ugr -- parps:")
#    print(parps)
    dots<-fnuser$dots # function does not know about fnuser ?? needed
    if (is.null(fnuser$gr)){
       stop("ugr: MUST specify a numerical derivative routine or character name")
    } else {
       usergr<-fnuser$gr
       if (is.character(usergr)) { # Character string means numerical gradient routine
          userfn<-fnuser$fn # Note that already have string. Must add function. 
          # ?? !! Note -- cannot use ufn because of scaling!
          # Consider fvalbest and parbest in opxfn (fnuser) structure and compare to see if usable.
          tgr<-try(tryg<-do.call(usergr, c(list(parps), list(userfn), dots)), silent=TRUE)
          fnuser$KFN<-npar+1+fnuser$KFN # Note we include the evaluation at base point
       } else {
          tgr<-try(tryg<-do.call("usergr", c(list(parps), dots)), silent = TRUE) 
          fnuser$KGR<-1+fnuser$KGR
       }
    }
    if ((class(tgr) == "try-error") || any(is.na(tryg)) || any(is.null(tryg)) || 
      any(is.infinite(tryg))) {
      tryg <- rep(.Machine$double.xmax, npar)
      attr(tryg, "inadmissible") <- TRUE
    }
    else {
      attr(tryg, "inadmissible") <- FALSE
    }
    if (any(is.null(tryg))) stop("NULL FUNCTION")
    if (fnuser$MAXIMIZE) tryg <- -tryg 
    tryg*fnuser$PARSCALE/fnuser$FNSCALE # handle the scaling
}
############# end ugr.R ##########################

