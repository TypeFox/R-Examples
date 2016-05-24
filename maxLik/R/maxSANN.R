maxSANN <- function(fn, grad=NULL, hess=NULL,
                    start, fixed = NULL,
                    control=NULL,
                    constraints = NULL,
                    finalHessian=TRUE,
                    parscale=rep(1, length=length(start)),
                    ... ) {
   ## Wrapper of optim-based 'SANN' optimization
   ## 
   ## contraints    constraints to be passed to 'constrOptim'
   ## finalHessian:   how (and if) to calculate the final Hessian:
   ##            FALSE   not calculate
   ##            TRUE    use analytic/numeric Hessian
   ##            bhhh/BHHH  use information equality approach
   ##
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values
   if(!inherits(control, "MaxControl")) {
      mControl <- maxControl(iterlim=10000L)
      mControl <- addControlList(mControl, control)
                           # default values
   }
   else {
      mControl <- control
   }
   mControl <- addControlList(mControl, list(...), check=FALSE)
   ## save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by 'optim( method="SANN" )')
   set.seed(slot(mControl, "sann_randomSeed"))
   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "SANN", fixed = fixed,
                      constraints = constraints,
                      finalHessian=finalHessian,
                      parscale = parscale,
                      control=mControl,
                      ... )
   return(result)
}

