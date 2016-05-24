"irfparF" <-
function (irfpar, lambdac, lambda, i, mudisp = FALSE, parmu=vector(), 
	 taudisp = FALSE, partau = vector(), dispmufun = "",
	 disptaufun = "", doublegaus = FALSE,multiplegaus=FALSE) 
{
	irfvec<-vector()
	if(mudisp){ 
		if(dispmufun == "discrete")
		        irfvec[1] <- parmu[i]
		#if(dispmufun == "logistic")
		#        irfvec[1] <- superLogistic(irfpar[1], parmu, lambda)
		if(dispmufun == "poly")
		        irfvec[1] <- simpPol(irfpar[1], parmu, lambda, lambdac)
	}
	else 
	     irfvec[1]<-irfpar[1]
	if(taudisp){
		if(disptaufun == "discrete")
			irfvec[2] <- partau[i]
		else {
		     if(dispmufun == "discrete")
			irfpar <- append(0, irfpar)
		     #if(disptaufun == "logistic")
		     #	  irfvec[2] <- superLogistic(irfpar[1], parmu, lambda)
		     if(disptaufun == "poly"){
			  irfvec[2] <- simpPol(irfpar[2], partau, lambda, 
				     lambdac)
		          
		     }
	        }
	}
	else {
	     irfvec[2] <- if(dispmufun != "discrete") irfpar[2]
			  else irfpar[1]
	}
	## assuming no dispersion for 2nd, 3rd, etc gaussians
	if(doublegaus) irfvec <- append(irfvec, irfpar[3:length(irfpar)])
  if(multiplegaus) irfvec <- append(irfvec, irfpar[3:length(irfpar)])
	irfvec
}

