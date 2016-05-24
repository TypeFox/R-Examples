vcov.segmented<-function (object, var.diff=FALSE, ...){
    if(inherits(object, "Arima")){
    v<-object$var.coef
    return(v)
    }
    if(inherits(object, "glm")){
        if(var.diff) warning("option var.diff=TRUE ignored with `glm' objects", call.=FALSE)
        so <- summary.glm(object, correlation = FALSE, ...)
        v<-so$dispersion * so$cov.unscaled
      } else {
        if(var.diff){
              if(length(object$nameUV$Z)>1) {
                var.diff<-FALSE
                warning("var.diff set to FALSE with multiple segmented variables", call.=FALSE)
                }
        v<-summary.segmented(object, var.diff=TRUE, correlation = FALSE, ...)$cov.var.diff
        } else {
          so<-summary.segmented(object, var.diff=FALSE, correlation = FALSE, ...)
          v<-so$sigma^2 * so$cov.unscaled
          }
        }
      return(v)
      }



