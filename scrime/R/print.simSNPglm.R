`print.simSNPglm` <-
function(x,...){
	beta<-format(c(x$beta0,x$beta))
	if(!is.null(x$prob)){
		sep<-"\n                   + "
		resp<-"logit(Prob(Y = 1)) ="
		spec<-if(is.null(x$p.cutoff)) "Draw from Bernoulli Distribution" 
			else paste("Prob(Y = 1) >",x$p.cutoff)	
	}
	else{
		sep<-"\n  + "
		resp<-"Y ="
		spec<-paste("Using",x$err.call)
	}
	model<-paste(beta[-1], " * (",x$ia,")",sep="",collapse=sep)	
	model<-paste(resp," ",beta[1],sep,model,sep="")
	if(!is.null(x$err))
		model<-paste(model,sep,if(x$beta0<0) " ","error",sep="")
	cat("Model:\n",model,"\n\n","Specification of Response: ",spec,"\n",sep="")
}

