mcmc.converge.check <-
function(model,factors,...){
	f1=c()
	for (f in factors){		
		pattern=paste('^\\w+:',f,'$',sep="")
		f1=append(f1,grep(pattern,names(posterior.mode(model$Sol))))
	}
	plot(model$Sol[,f1])
}
