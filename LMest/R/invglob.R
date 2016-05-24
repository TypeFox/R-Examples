invglob <-
function(eta){

	b = length(eta)
	if(b == 1) p = c(1,exp(eta))/(1+exp(eta))
	else{
		pc = c(1/(1+exp(eta)),1)
	  	Di = cbind(rep(0,b),diag(b))-cbind(diag(b),rep(0,b))
	  	Di = rbind(c(1,rep(0,b)), Di)
	  	p = Di%*%pc
	}
}
