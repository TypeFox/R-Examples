bread.mlm <-
function(x, ...) {
   d <- length(coef(x))
   rval <- diag(d)
   colnames(rval) <- rownames(rval) <- colnames(coef(x))
   return(rval)
 }

