compare.etasclass <-
function(etas1,etas2){
# compare two etasclass objects
#check for classes
	  if(class(etas1)!="etasclass")stop("first argument must be an etasclass object")
	  if(class(etas2)!="etasclass")stop("second argument must be an etasclass object")
	  ### check for comparability of the two input objects
	  ## same threshold
	  ## same catalog
	  ## same domain (space time)
	  ##
	  params	=0
	  npar.est	=etas1$params.ind+etas2$params.ind
	  sqm		=sqrt((etas1$sqm^2+etas2$sqm^2)/npar.est)
	  if(sqm>0) params=(etas1$params-etas2$params)/sqm
	  AIC=min(etas1$AIC)-min(etas2$AIC)
	  if(etas1$onlytime&etas2$onlytime){
	  weights	=0
	  weights.std	=0
	  cor.weights	=0
	  }
	  {
	  weights	=etas1$rho.weights-etas2$rho.weights
	  weights.std	=(etas1$rho.weights-etas2$rho.weights)/sqrt((etas1$rho.weights*(1-etas1$rho.weights)+etas2$rho.weights*(1-etas2$rho.weights))/2)
	  cor.weights	=cor(etas1$rho.weights,etas2$rho.weights)
	  }
	  cor.trig	=cor(etas1$l-etas1$params[1]*etas1$back.dens,etas2$l-etas2$params[1]*etas2$back.dens)
	  cor.back	=cor(etas1$back.dens,etas2$back.dens)
return(list(params=params,AIC=AIC,weights=weights,weights.std=weights.std,cor.weights=cor.weights,cor.trig=cor.trig,cor.back=cor.back))
  
	  }
