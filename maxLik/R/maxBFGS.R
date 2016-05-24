maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start, fixed = NULL,
                    control=NULL,
                    constraints=NULL,
                    finalHessian=TRUE,
                    parscale=rep(1, length=length(start)),
                    ## sumt parameters
                    ...) {
   ## Wrapper of optim-based 'BFGS' optimization
   ## 
   ## contraints    constraints to be passed to 'constrOptim'
   ## finalHessian:   how (and if) to calculate the final Hessian:
   ##            FALSE   not calculate
   ##            TRUE    use analytic/numeric Hessian
   ##            bhhh/BHHH  use information equality approach
   ##
   ## ...           further arguments to fn() and grad()
   if(!inherits(control, "MaxControl")) {
      mControl <- addControlList(maxControl(iterlim=200), control)
                           # default values
   }
   else {
      mControl <- control
   }
   mControl <- addControlList(mControl, list(...), check=FALSE)
   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "BFGS", fixed = fixed,
                      constraints = constraints,
                      finalHessian=finalHessian,
                      parscale = parscale,
                      control=mControl,
      ... )

   return(result)
}
