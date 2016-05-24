inv_glob <-function(eta,type="g",der=F){

# if only one logit
	D = NULL
	eta = as.vector(eta)
  if(length(eta)==1){
  	p = exp(eta)
  	p = c(1,p)/(1+p)
  	out = list(p=p,D=D)
  }else{
# Invert global logits collected in column vectoreta
	  if(length(eta)==1){
  		p = exp(eta); p = p/(1+p)
  		p = c(1-p,p)
	   	out = list(p=p,D=D)
	  }else{
# invert parametrization
		if(type=="g"){
			eeta = exp(eta)
			pc = c(0,1/(1+eeta),1)
			p = diff(pc)
			if(der){
				zz = rep(0,length(eta))
				D = diff(rbind(zz,diag(-eeta/(1+eeta)^2),zz))
			}
		}
		if(type=="l"){
			pc = exp(c(0,cumsum(eta)))
    		p = pc/sum(pc)
    	}
    	out = list(p=p,D=D)
    }
    }
}
