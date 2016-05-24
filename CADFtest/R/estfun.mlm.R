estfun.mlm <-
function(x, ...) {
   psi <- t(t(model.response(model.frame(x))) - as.vector(coef(x)))
   colnames(psi) <- colnames(coef(x))
   return(psi)
 }

