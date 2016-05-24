`coef.eRm` <-
function(object, parm = "beta", ...) {         # option "beta" added rh 2010-03-07
   if(parm == "beta")
       object$betapar
   else if(parm == "eta")
       object$etapar
   else
       stop("'parm' incorrectly specified")
}
