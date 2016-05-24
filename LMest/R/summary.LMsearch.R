summary.LMsearch<-function(object,...){ 
	cat("Call:\n")
	print(object$call)
	kdim = length(object$kv)
	np = lk = AIC = BIC= rep(0,kdim)
	cont=0
	for(k in object$kv){
		cont=cont+1
		np[cont] = object$out.single[[k]]$np
		lk[cont] = object$lkv[k]
		AIC[cont] = object$aicv[k]
		BIC[cont] = object$bicv[k]		
	}
	print(cbind(states=object$kv,lk=lk,np=np,AIC=AIC,BIC=BIC))
	
}