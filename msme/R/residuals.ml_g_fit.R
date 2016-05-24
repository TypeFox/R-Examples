residuals.ml_g_fit <-
   function(object, type = c("raw","ss"), ...) {
   type <- match.arg(type)
   e.hat <- object$y - fitted(object)
   if (type == "ss") {
      e.hat <- e.hat /
              (object$sigma.hat * sqrt(1 - diag(hatvalues(object)))) 
    }
   return(e.hat)
}
