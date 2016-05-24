#' @importFrom numDeriv grad
optimx.setup <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {
### To return in optcfg: fname, npar ??, method, ufn, ugr, ctrl, have.bounds

# Get real name of function to be minimized
  fname<-deparse(substitute(fn))
  if (!is.null(control$trace) && control$trace>0) {
	cat("fn is ",fname,"\n")
  }
  optcfg<-list()
  optcfg$fname<-fname

# Only one ref to parameters -- to get npar here
  npar <- length(par) # !! NOT CHECKED in case par not well-defined
  optcfg$npar <- npar

# Set control defaults
    ctrl <- list(
	follow.on=FALSE, 
	save.failures=TRUE,
	trace=0,
	kkt=TRUE,
	all.methods=FALSE,
	starttests=TRUE,
	maximize=FALSE,
	dowarn=TRUE, 
        usenumDeriv=FALSE,
	kkttol=0.001,
	kkt2tol=1.0E-6,
	badval=(0.5)*.Machine$double.xmax,
	scaletol=3
    ) 
    
# Note that we do NOT want to check on the names, because we may introduce 
#    new names in the control lists of added methods
#    if (!all(namc %in% names(ctrl))) 
#        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
# However, we do want to substitute the appropriate information. 
# removed copy of hessian to control$kkt
    ncontrol <- names(control)
    nctrl <- names(ctrl)
    for (onename in ncontrol) {
       if (onename %in% nctrl) {
           ctrl[onename]<-control[onename]
       } else {
           ctrl[onename]<-control[onename]
       }
    }
    if (is.null(control$kkt)) { # turn off kkt for large matrices
      ctrl$kkt<-TRUE # default it to compute KKT tests
      if (is.null(gr)) { # no analytic gradient
         if (npar > 50) {
           ctrl$kkt=FALSE # too much work when large number of parameters
           if (ctrl$trace>0) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      } else {
         if (npar > 500) {
            ctrl$kkt=FALSE # too much work when large number of parameters, even with analytic gradient
            if (ctrl$trace>0) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      }
    } else { # kkt is set
      if (control$kkt) {
        if (is.null(gr)) {
           if (npar > 50) {
             if ((ctrl$trace>0) && ctrl$dowarn) warning("Computing hessian for gr NULL, npar > 50, can be slow\n")
           }
        } else {
           if (npar > 500) {
             if ((ctrl$trace>0) && ctrl$dowarn) warning("Computing hessian with gr code, npar > 500, can be slow\n")
           }
        }
      }
    }
    optcfg$ctrl <- ctrl
# reset the function if we are maximizing
  ufn <- fn
  ugr <- gr
  uhess <- hess
  if ((! is.null(control$maximize)) && control$maximize ) { 
        cat("Maximizing -- use negfn and neggr\n")
        if (! is.null(control$fnscale)) { 
 		stop("Mixing controls maximize and fnscale is dangerous. Please correct.")
        } # moved up 091216
        optcfg$ctrl$maximize<-TRUE
        ufn <- function (par, ...) { # negate the function for maximizing
	   val<-(-1.)*fn(par,...)
        } # end of ufn = negfn
        if (! is.null(gr)) { 
           ugr <- function(par, userfn=ufn, ...) {
               gg <- (-1)*gr(par, ...)
           }
        } else { ugr <- NULL } # ensure it is defined
        if (! is.null(hess) ) {
           uhess <- function(par, ...) {
               hh <- (-1)*hess(par, ...)
           }
        } else { uhess <- NULL } # ensure it is defined
  } else { 
     optcfg$ctrl$maximize <- FALSE # ensure defined
  } # define maximize if NULL
  optcfg$usenumDeriv<-FALSE # JN130703
  if (is.null(gr) && ctrl$dowarn && ctrl$usenumDeriv) {
     warning("Replacing NULL gr with 'numDeriv' approximation")
     optcfg$usenumDeriv<-TRUE
     ugr <- function(par, userfn=ufn, ...) { # using grad from numDeriv
        tryg<-numDeriv::grad(userfn, par, ...)
     } # Already have negation in ufn if maximizing
  }
  optcfg$ufn <- ufn
  optcfg$ugr <- ugr
  optcfg$uhess <- uhess

# Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) { have.bounds<-TRUE # set this for convenience
  } else { have.bounds <- FALSE }
  optcfg$have.bounds <- have.bounds

  # Check that we have the functions we need
#   if (! require(numDeriv, quietly=TRUE) ) stop("Install package `numDeriv'", call.=FALSE)

  # List of methods in base or stats, namely those in optim(), nlm(), nlminb()
  bmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb")
# SANN has no termination for optimality, only a maxit count for
#    the maximum number of function evaluations; remove DEoptim for now -- not useful 
#    for smooth functions. Code left in for those who may need it.
  # List of methods in packages. 
  allmeth <- bmeth
# Now make sure methods loaded
   allmeth <- bmeth # start with base methods
   testload <- TRUE # This is a temporary fix for NAMESPACE changes in R 3.1.2
#   testload <- suppressWarnings(require(BB, quietly=TRUE))
   if (testload)  allmeth<-c(allmeth,"spg")
   else if (ctrl$trace>0) { warning("Package `BB' Not installed", call.=FALSE) }

#   testload <- suppressWarnings(require(ucminf, quietly=TRUE))
   if (testload)  allmeth<-c(allmeth,"ucminf")
   else if (ctrl$trace>0) { warning("Package `ucminf' Not installed", call.=FALSE) }
   
#   testload <- suppressWarnings(require(Rcgmin, quietly=TRUE))
   if (testload)  allmeth<-c(allmeth,"Rcgmin")
   else if (ctrl$trace>0) { warning("Package `Rcgmin' Not installed", call.=FALSE) }
   
#   testload <- suppressWarnings(require(Rvmmin, quietly=TRUE))
   if (testload)  allmeth<-c(allmeth,"Rvmmin")
   else if (ctrl$trace>0) { warning("Package `Rvmmin' Not installed", call.=FALSE) }
   
#   testload <- suppressWarnings(require(minqa, quietly=TRUE))
   if (testload) { allmeth<-c(allmeth, "newuoa", "bobyqa")  }
   else if (ctrl$trace>0) { warning("Package `minqa' (for uobyqa, newuoa, and bobyqa) Not installed", call.=FALSE) }
   # leave out uobyqa in CRAN version 120421 (from earlier 1104 change)

#   testload <- suppressWarnings(require(dfoptim, quietly=TRUE))
   if (testload)  allmeth<-c(allmeth,"nmkb", "hjkb")
   else if (ctrl$trace>0) { warning("Package `dfoptim' Not installed", call.=FALSE) }
   
   bdsmeth<-c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "bobyqa", "nmkb", "hjkb")
  # Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) allmeth <- allmeth[which(allmeth %in% bdsmeth)]
  if (ctrl$all.methods) { # Changes method vector!
	method<-allmeth
        if (ctrl$trace>0) {
		cat("all.methods is TRUE -- Using all available methods\n")
		print(method)
	}
  } 

  # Partial matching of method string allowed
  # avoid duplicates here
  # 2011-1-17 JN: to set L-BFGS-B
  method <- try(unique(match.arg(method, allmeth, several.ok=TRUE) ),silent=TRUE)
  if (class(method)=="try-error") {
     warning("optimx: No match to available methods")
     method<-NULL
     nmeth<-0
  } else {
     nmeth <- length(method) # number of methods requested
  } # JN 2011-1-17 fix for default when there are bounds
  if ((nmeth==0) && have.bounds) {
      method="L-BFGS-B"
      if (ctrl$dowarn) warning("Default method when bounds specified is L-BFGS-B to match optim()")
      nmeth<-1
  }
  ## Check that methods are indeed available and loaded
  for (i in 1:nmeth) {
     cmeth <- method[i]
     if (ctrl$trace > 0) cat("Looking for method = ",cmeth,"\n")
     if (! (cmeth %in% allmeth) ) {
         errmsg <- paste(cmeth," not found in any list of methods available")
         stop(errmsg, call.=FALSE)
     } # otherwise the method is available, and just needs to be loaded
  } # end check methods available
  if (ctrl$trace>1) {
    cat("Methods to be used:")
    print(method)
  }
  optcfg$method <- method
  optcfg # return the structure
} ## end of optimx.setup


