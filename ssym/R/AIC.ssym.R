AIC.ssym <-
function(object, ...){
	gle <- 0
	if(object$p>0) gle <- gle + object$p
	if(sum(object$qm) > 0)	gle <- gle + sum(object$gle.mu)
	if(object$l>0) gle <- gle + object$l
	if(sum(object$q) > 0) gle <- gle + sum(object$gle.phi)

    AIC <- round(-2*sum(object$lpdf) + 2*(gle), digits=3)
    y <- object$z_es*sqrt(object$phi.fitted) + object$mu.fitted
	if(object$censored==FALSE) attr(AIC,"log") <- round(-2*sum(object$lpdf) + 2*(gle) + 2*sum(y), digits=3)
	else attr(AIC,"log") <- round(-2*sum(object$lpdf) + 2*(gle) + 2*sum(y[object$event==0]), digits=3)
    AIC
}
