residuals.ssym <-
function(object, ...){
if(object$censored=="FALSE")
list(mu=sqrt(object$deviance.mu)*ifelse(object$z_es>=0,1,-1), phi=sqrt(object$deviance.phi)*ifelse(object$z_es>=0,1,-1), overall=qnorm(object$cdfz), ordinary=object$z_es)
else
list(mu=sqrt(object$deviance.mu)*ifelse(object$z_es>=0,1,-1), phi=sqrt(object$deviance.phi)*ifelse(object$z_es>=0,1,-1))}
