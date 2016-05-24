#' @importFrom stats glm ppois dpois poisson
poisson.scores <- function(y, X, ...){
  N <-  length(y)
  mod <- glm(y~X, family=poisson(), ...)

  dl.dtheta <- cbind(y- mod$fitted.value, (y- mod$fitted.value)*X)
  d2l.dtheta.dtheta <- - crossprod(cbind(1, X)*sqrt(mod$fitted))
  
  presid <- ppois(y-1, mod$fitted.values) + ppois(y, mod$fitted.values) -1
  
  #dpresid.dlambda <- 2*(ppois(y-2, mod$fitted.values) - ppois(y-1, mod$fitted.values)) + dpois(y-1, mod$fitted.values) - dpois(y, mod$fitted.values)
  #### can be further simplified as
  dpresid.dlambda <- -  (dpois(y-1, mod$fitted.values) + dpois(y, mod$fitted.values))
  dlambda.dtheta <- mod$fitted.values * cbind(1, X)
  dpresid.dtheta <- t(dpresid.dlambda * dlambda.dtheta)
  
  pearson.resid <- (y - mod$fitted.values)/sqrt(mod$fitted.values)
  dpearson.resid.dlambda <- - 0.5 * ( mod$fitted.values^(-1/2) + mod$fitted.values^(-3/2)*y)
  dpearson.resid.dtheta <- t(dpearson.resid.dlambda * dlambda.dtheta)
  
  list(mod = mod,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       presid = presid,
       dpresid.dtheta = dpresid.dtheta,
       pearson.resid = pearson.resid,
       dpearson.resid.dtheta = dpearson.resid.dtheta)
}
