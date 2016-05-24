summary.GORH <-
function(object,...){
      coef<-object$ParEst$Beta
	se<-sqrt(diag(object$ParVcov$vcov.bg))
	tval<-coef/se
	TAB<-cbind(Estimate=coef,
           StdErr=se,
           t.value=tval,
           p.value=2*pnorm(-abs(tval)))
	rownames(TAB)<-all.vars(object$survfun)[-c(1:2)]
	res<-list(call=object$call,
          coefficients=TAB,loglik=round(object$ParEst$loglik,2))
	class(res)<-"summary.GORH"
	return(res)
}
