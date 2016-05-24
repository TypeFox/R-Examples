optimx.check <- function(par, ufn, ugr, uhess, lower=-Inf, upper=Inf, 
             hessian=FALSE, ctrl, have.bounds=FALSE, usenumDeriv=FALSE, ...) {
##            method=NULL, itnmax=NULL, hessian=FALSE,
##            ctrl=list(),...) {

## Should be run whenever we are not sure parameters and function are
## admissible. 
##
## Inputs: ?? 
##   control list -- ctrl. ufn, ugr

## Outputs: ?? failed-checks info.

###############################################################################

## Code more or less common to funtest, funcheck and optimx <<<
# Check parameters are in right form
  if (!is.null(dim(par))) stop("Parameter should be a vector, not a matrix!", call. = FALSE)
  if (! is.vector(par) ) {
	stop("The parameters are NOT in a vector")
  }
  npar<-length(par)
  optchk<-list() # output of the checks
  if (ctrl$starttests) {
	# Check parameters in bounds (090601: As yet not dealing with masks ??)
	#  bdmsk<-as.vector(bdmset[k, ])
	infeasible<-FALSE
	if (ctrl$trace > 0) cat("Function has ",npar," arguments\n")
	if (have.bounds) {
    	  # Expand bounds to vectors if needed
          # Note 20100610: we do not check if there is a vector of wrong length.??
    	  if (length(lower)==1 ) lower <- rep(lower, npar)
    	  if (length(upper)==1 ) upper <- rep(upper, npar)
    	  bstate<-vector(mode="character", length=npar)
    	  for (i in 1:npar) {
     	    if ( (lower[i]<=par[i]) && (par[i]<=upper[i])) {
    	      bstate[i]<-" In Bounds "
            } else { 
            #   if (bdmsk[i]!=0) {
                  infeasible<-TRUE
            #   }
                if (lower[i]>par[i]) {bstate[i]<-" Out of Bounds LOW" } else { bstate[i]<-" Out of Bounds HIGH " }
            } # end if in bounds
#           if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bdmsk[i],"   ",bstate,"\n")
            if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bstate,"\n")
          } # end of for loop over parameter vector elements
	  if (infeasible) { ## ?? maybe don't want to stop ??
        	stop("Infeasible point, no further tests")
	  } 
  	} # end have.bounds
        # Check if function can be computed
        firsttry<-try(finit<-ufn(par, ...), silent=TRUE ) # 20100711
        # Note: This incurs one EXTRA function evaluation because optimx is a wrapper for other methods
        if (class(firsttry) == "try-error") {
    	   infeasible <- TRUE
           stop("Cannot evaluate function at initial parameters")
        }
        # Also check that it is returned as a scalar
       if (!(is.vector(finit) && (length(finit)==1)) || is.list(finit) || 
            is.matrix(finit) || is.array(finit) || ! is.numeric(finit) ) {
           stop("Function provided is not returning a scalar number")
       }
       if (is.infinite(finit) || is.na(finit)) {
          stop("Function returned is infinite or NA (non-computable)")
       }
  }


  if (ctrl$starttests) {
     optchk$grbad <- FALSE
     if (! is.null(ugr) && ! usenumDeriv){ # check gradient
       gname <- deparse(substitute(ugr))
       if (ctrl$trace>0) cat("Analytic gradient from function ",gname,"\n\n")
          fval <- ufn(par,...) 
          gn <- grad(func=ufn, x=par,...) # 
          ga <- ugr(par, ...)
          # Now test for equality (090612: ?? There may be better choices for the tolerances.
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(gn-ga))/(1 + abs(fval)) >= teps) {
            stop("Gradient function might be wrong - check it! \n", call.=FALSE)
            optchk$grbad <- TRUE # Never get here if we stop ??
          }
       } else if (ctrl$trace>0) cat("Analytic gradient not made available.\n")

       optchk$hessbad <- FALSE
       if (! is.null(uhess)){ # check Hessian
          hname <- deparse(substitute(uhess))
          if (ctrl$trace>0) cat("Analytic hessian from function ",hname,"\n\n")
          hn <- hessian(func=ufn, x=par,...) # ?? should we use dotdat
          ha <- uhess(par, ...)
          # Now test for equality
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) stop("Hessian function might be wrong - check it! \n", call.=FALSE)
          optchk$hessbad <- TRUE
       } else if (ctrl$trace>0) cat("Analytic Hessian not made available.\n")
   }
# Scaling check  091219
    if (ctrl$starttests) {
        optchk$scalebad <- FALSE
	srat<-scalecheck(par, lower, upper,ctrl$dowarn)
	sratv<-c(srat$lpratio, srat$lbratio)
	if (max(sratv,na.rm=TRUE) > ctrl$scaletol) { 
		warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
		if (ctrl$dowarn) warning(warnstr)
             optchk$scalebad <- TRUE
	}
        if (ctrl$trace>0) {
		cat("Scale check -- log parameter ratio=",srat$lpratio,"  log bounds ratio=",srat$lbratio,"\n")
	}
    }
# end scaling check

## ?? what to return
    optchk
} ## end of optimx.check

