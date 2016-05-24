`spatial.location` <-
function(X, score=c("sign","signrank"), init=NULL, shape=TRUE, steps=Inf, maxiter=500, eps=1e-6, na.action=na.fail)
{ 
 X <- na.action(X)

 score<-match.arg(score)
 switch(score,
       "sign"=
       ae.spatial.median(X, init=init, shape=shape, steps=steps, 
       maxiter=maxiter, eps=eps) 
       ,
       "signrank"=
       ae.hl.estimate(X, init=init, shape=shape, steps=steps,
       maxiter=maxiter, eps=eps)
       )
}

