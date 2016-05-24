summary.TransModel <-
function(object,...){
	se<-sqrt(diag(object$vcov))
	tval<-coef(object)/se
	TAB<-cbind(Estimate=coef(object),
           StdErr=se,
           t.value=tval,
           p.value=2*pt(-abs(tval),df=object$df))
      col<-ncol(object$data)-1
	rownames(TAB)<-all.vars(object$formula)[-c(1:2)]
	res<-list(call=object$call,
          coefficients=TAB)
	class(res)<-"summary.TransModel"
	return(res)
}
