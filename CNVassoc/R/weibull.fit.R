weibull.fit <-
function(x, y, cens, weights){
  if (any(weights <= 0)){
    keep <- weights > 0
    x <- x[keep, ]
    y <- y[keep]
    cens <- cens[keep] 
    weights <- weights[keep]    
  }
  fit<-survreg(Surv(y,cens)~x-1, weights=weights)
  gamma<-coef(fit)
  alpha<-1/fit$scale
  beta<- -gamma*alpha
  return(list(beta=beta,alpha=alpha,mm=model.matrix(fit)))
}