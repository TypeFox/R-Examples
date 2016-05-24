
## similar to residuals.gam()....

residuals.scam <-function(object, type = c("deviance", "pearson","scaled.pearson", "working", "response"),...)
# calculates residuals for scam object the same as residulas.gam()...
{ type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  y.mu <- y-mu
##  family <- object$family
  wts <- object$prior.weights
  res<- switch(type,working = object$residuals,
         scaled.pearson = y.mu*wts^.5/sqrt(object$sig2*object$family$variance(mu)),
              pearson = y.mu*wts^.5/sqrt(object$family$variance(mu)),
              deviance = { d.res <- sqrt(pmax(object$family$dev.resids(y,mu,wts),0))
                           ifelse(y>mu , d.res, -d.res)             
                         },
              response = y.mu)
  res <- naresid(object$na.action,res)
  res
}

