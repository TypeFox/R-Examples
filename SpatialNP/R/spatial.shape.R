`spatial.shape` <-
function(X, score=c("sign","symmsign","rank","signrank"), fixed.loc=FALSE, location=NULL, init=NULL, steps=Inf, eps=1e-06, maxiter=100, na.action=na.fail)
{ 
 X <- na.action(X)

 score<-match.arg(score)
 switch(score,
  "sign" = signs.shape(X, fixed.loc=fixed.loc, location=location, init=init, steps=steps, eps=eps, maxiter=maxiter),

  "symmsign" = symmsign.shape(X, init=init, steps=steps, eps=eps, maxiter=maxiter),

  "rank" = rank.shape(X, init=init, steps=steps, eps=eps, maxiter=maxiter),
 
  "signrank" = signrank.shape(X, fixed.loc=fixed.loc, location=location, init=init, steps=steps, eps=eps, maxiter=maxiter),
 )
}

