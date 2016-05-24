maxCG <- function(fn, grad=NULL, hess=NULL,
                  start, fixed = NULL,
                  control=NULL,
                  constraints=NULL,
                  finalHessian=TRUE,
                  parscale=rep(1, length=length(start)),
                  ...) {
   ## Wrapper of optim-based 'Conjugate Gradient' optimization
   ## 
   ## contraints    constraints to be passed to 'constrOptim'
   ## hessian:   how (and if) to calculate the final Hessian:
   ##            FALSE   not calculate
   ##            TRUE    use analytic/numeric Hessian
   ##            bhhh/BHHH  use information equality approach
   ## ... :      further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values
   ##
   if(!inherits(control, "MaxControl")) {
      mControl <- addControlList(maxControl(iterlim=500), control)
                           # default values
   }
   else {
      mControl <- control
   }
                           # default, user values
   mControl <- addControlList(mControl, list(...), check=FALSE)
                           # open values
   result <- maxOptim( fn = fn, grad = grad, hess = hess,
                      start = start, method = "CG", fixed = fixed,
                      constraints = constraints,
                      finalHessian=finalHessian,
                      parscale = parscale,
                      control=mControl,
                      ... )

   return(result)
}
