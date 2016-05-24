vcov.ssym <-
function(object, ...){
if(object$censored==FALSE) list(mu=object$vcov.mu, phi=object$vcov.phi)
else list(mu=object$vcov.mu, phi=object$vcov.phi, cov.mu.phi=object$cov.mu.phi)}
