#' Fit detection function using key-adjustment functions
#'
#' Fit detection function to observed distances using the key-adjustment
#' function approach.  If adjustment functions are included it will alternate
#' between fitting parameters of key and adjustment functions and then all
#' parameters much like the approach in the CDS and MCDS Distance FORTRAN code.
#' To do so it calls \code{detfct.fit.opt} which uses the R optim function
#' which does not allow non-linear constraints so inclusion of adjustments does
#' allow the detection function to be non-monotone.
#'
#' @aliases detfct.fit
#' @param ddfobj detection function object
#' @param optim.options control options for optim
#' @param bounds bounds for the parameters
#' @param misc.options miscellaneous options
#' @return fitted detection function model object with the following list
#'   structure \item{par}{final parameter vector} \item{value}{final negative
#'   log likelihood value} \item{counts}{number of function evaluations}
#'   \item{convergence}{see codes in optim} \item{message}{string about
#'   convergence} \item{hessian}{hessian evaluated at final parameter values}
#'   \item{aux}{ a list with 20 elements \itemize{ \item maxit: maximum number
#'   of iterations allowed for optimization \item lower: lower bound values for
#'   parameters \item upper: upper bound values for parameters \item setlower:
#'   TRUE if they are user set bounds \item setupper: TRUE if they are user set
#'   bounds \item point: TRUE if point counts and FALSE if line transect \item
#'   int.range: integration range values \item showit: integer value that
#'   determines information printed during iteration \item integral.numeric
#'   if TRUE compute logistic integrals numerically \item
#'   breaks: breaks in distance for defined fixed bins for analysis \item
#'   maxiter: maximum iterations used \item refit: if TRUE, detection function
#'   will be fitted more than once if parameters are at a boundary or when
#'   convergence is not achieved \item nrefits: number of refittings \item
#'   \item mono: if TRUE montonicity will be enforced \item
#'   mono.strict: if TRUE, then strict monotonicity is enforced; otherwise weak
#'   \item width: radius of point count or half-width of strip \item
#'   standardize: if TRUE, detection function is scaled so g(0)=1 \item ddfobj:
#'   distance detection function object; see \code{\link{create.ddfobj}} \item
#'   bounded: TRUE if parameters ended up a boundary (I think) \item model:
#'   list of formulas for detection function model (probably can remove this)
#'   }}
#' @author Dave Miller; Jeff Laake
detfct.fit <- function(ddfobj,optim.options,bounds,misc.options){
  # Functions Used: assign.par, detfct.fit.opt, get.par

  # show debug information
  showit <- misc.options$showit

  # How small is small?
  epsilon <- sqrt(.Machine$double.eps)

  # keep a history of how the optimisation is doing
  # stores: convergence status (0=GOOD), lnl, pars
  misc.options$optim.history <- rep(NA, length(getpar(ddfobj))+2)

  # Count how we're doing...
  iter <- 0
  metaiter <- 0

  # If we have no adjustments then we can just do the optimisation -- no
  # complicated stuff. Also if:
  #  - we have uniform detection function
  #  - we're enforcing monotonicity
  if(is.null(ddfobj$adjustment) | ddfobj$type=="unif" |
     misc.options$mono | misc.options$nofit){

    if(misc.options$mono & ddfobj$type!="unif"){
      ## get best key pars first, not enforcing monotonicity
      save.mono <- misc.options$mono
      save.mono.strict <- misc.options$mono.strict
      misc.options$mono <- FALSE
      misc.options$mono.strict <- FALSE
      lt <- detfct.fit.opt(ddfobj, optim.options, bounds,
                           misc.options, fitting="key")
      misc.options$mono <- save.mono
      misc.options$mono.strict <- save.mono.strict
      # assign those parameters
      ddfobj <- assign.par(ddfobj, lt$par)
      # recalculate lower/upper bounds
      if(check.bounds(lt, bounds$lower, bounds$upper, ddfobj,
                            showit,bounds$setlower, bounds$setupper)){
        bounds <- setbounds(rep(NA, length(lt$par)), rep(NA, length(lt$par)),
                            lt$par, ddfobj)
      }
    }

    lt <- detfct.fit.opt(ddfobj, optim.options, bounds, misc.options)

  }else{
  # Otherwise we need to play around...

    # think this needs to live elsewhere, but let's leave it here for
    # the moment
    if(!is.null(ddfobj$adjustment) && ddfobj$adjustment$series=="herm")
      ddfobj$adjustment$parameters <- rep(1, length(ddfobj$adjustment$order))

    initialvalues <- getpar(ddfobj)
    # This holds the previous values, to test for convergence
    lastvalues <- initialvalues

    # Now start to alternate between adjustment and key
    if(showit>=2)
      cat("DEBUG: starting to cycle the optimisation...\n")

    # Fudge to make this work the first time.
    firstrun<-TRUE

    while((iter < misc.options$maxiter) &&
           (all((abs(initialvalues-lastvalues)/(epsilon+abs(lastvalues))) <
               (epsilon/sqrt(nrow(ddfobj$xmat)))) |
              firstrun)){

      # Variable to count sub iterations... :)
      metaiter<-0
      firstrun<-FALSE

      # loop through fitting the adjustment, key and full detection function
      for(fitting in c("adjust","key","all")){

        # don't do refitting when we are just fitting key or adjustments
        if(fitting == "adjust" | fitting=="key"){
          refit.save<-misc.options$refit
          misc.options$refit<-FALSE
        }

        if(showit >= 2) {
          cat("DEBUG:",fitting,"iteration ",iter,".",metaiter,"\n")
          cat("DEBUG: initial values =",
               paste(round(initialvalues, 7), collapse=", "),"\n")
        }

        lt <- try(detfct.fit.opt(ddfobj,optim.options,bounds,misc.options,
                             fitting=fitting))
        metaiter <- metaiter+1

        # report failure
        if(all(class(lt)=="try-error")){
          if(showit>=2){
            cat("DEBUG: iteration",iter,", fitting", fitting,"failed.\n")
          }
          if(fitting=="all"){
            stop("Fitting failed! Try again with better initial values.\n")
          }
        }else{
          # update bounds
          bounds<-lt$bounds

          if(showit >= 2){
            cat("DEBUG: iteration", paste(iter,metaiter,sep="."),
                ":\n       Converge   =",lt$converge,
                "\n       lnl        =",lt$value,
                "\n       parameters =",paste(round(lt$par,7),collapse=", "),"\n")
          }

          # were any of the pars NA?
          # if so, reset to initialvalues
          if(any(is.na(lt$par))){
            lt$par[is.na(lt$par)]<-initialvalues[is.na(lt$par)]
          }
          # Rebuild initialvalues again
          ddfobj <- assign.par(ddfobj,lt$par)
          initialvalues <- getpar(ddfobj)

          # save the optimisation history
          optim.history <- lt$optim.history
          misc.options$optim.history <- optim.history
        }

        # restore the refit status
        misc.options$refit<-refit.save

      }

      iter<-iter+1
    }

    if(iter>misc.options$maxiter){
      warning("Maximum iterations exceeded!",
          iter,">",misc.options$maxiter,"\n")
    }

  }

  if(showit>=1){
    cat("\nDEBUG: Convergence!",
        "\n       Iteration ",paste(iter,metaiter,sep="."),
        "\n       Converge   =",lt$converge,
        "\n       lnl        =",lt$value,
        "\n       parameters =",paste(round(lt$par,7),collapse=", "),"\n")
  }

  # get rid of the first (dummy) line of the optimisation history
  lt$optim.history<-lt$optim.history[-1,]

  # Return some (hopefully correct) results
  return(lt)
}
