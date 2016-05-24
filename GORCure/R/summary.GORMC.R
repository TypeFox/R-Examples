summary.GORMC <-
function(object,...){
      coef<-c(object$ParEst$Eta,object$ParEst$Beta)
	se<-sqrt(diag(object$ParVcov$vcov.bg))
	tval<-coef/se
	TAB<-cbind(Estimate=coef,
           StdErr=se,
           t.value=tval,
           p.value=2*pnorm(-abs(tval)))
	rownames(TAB)<-c(paste("INC:",c("Intercept",all.vars(object$curefun))),paste("LAT:",all.vars(object$survfun)[-c(1:2)]))
	res<-list(call=object$call,
          coefficients=TAB,loglik=round(object$ParEst$loglik,2))
	class(res)<-"summary.GORMC"
	return(res)
}
