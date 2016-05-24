logLik.ssym <-
function(object, ...){
logLik <- round(sum(object$lpdf),digits=3)
    y <- object$z_es*sqrt(object$phi.fitted) + object$mu.fitted
	if(object$censored==FALSE) attr(logLik,"log") <- round(sum(object$lpdf) - sum(y), digits=3)
	else attr(logLik,"log") <- round(sum(object$lpdf) - sum(y[object$event==0]), digits=3)
    logLik
}
