UpdateRP.FLrw <-
function(survObj, priorPara, mcmcPara, ini){
	
	n 				<- survObj$n
	p				<- survObj$p
	x				<- survObj$x	
	
	J				<- priorPara$J
	ind.r			<- priorPara$ind.r
	ind.d			<- priorPara$ind.d
	ind.r_d			<- priorPara$ind.r_d		
	numBeta			<- mcmcPara$numBeta	
	beta.prop.me	<- mcmcPara$beta.prop.me
	beta.prop.var	<- mcmcPara$beta.prop.var
	
	xbeta			<- ini$xbeta
	be.ini			<- ini$beta.ini
	h				<- ini$h
	sigmaSq			<- ini$sigmaSq
	tauSq			<- ini$tauSq
	wSq				<- ini$wSq
	#sd.be			<- ini$sd.be
	
	
	updatej	<- sample(1:p, numBeta)
	accept 	<- rep(0, p)
	
			
	for(j in updatej){
		
		be.prop <- be.ini

		xbeta[xbeta > 700] <- 700
		exp.xbeta		<- exp(xbeta)
		exp.xbeta.mat	<- matrix(rep(exp.xbeta, J), n, J)
		
		first.sum 	<- colSums(exp.xbeta.mat*ind.r_d)

		h.mat					<- matrix(rep(h, n), n, J, byrow = T)
		h.exp.xbeta.mat 	<- - h.mat * exp.xbeta.mat
		h.exp.xbeta.mat[h.exp.xbeta.mat > -10^(-7)] <- -10^(-7)
		second.sum 	<- colSums(log(1 - exp(h.exp.xbeta.mat))*ind.d)

		loglh.ini <- sum(-h*first.sum + second.sum)
		
		be.prop[j] <- rnorm(1, mean = beta.prop.me[j], sd = sqrt(beta.prop.var))

		xbeta.prop			<- xbeta - x[,j]*be.ini[j] + x[,j]*be.prop[j] 
		xbeta.prop[xbeta.prop > 700] <- 700
		exp.xbeta.prop		<- exp(xbeta.prop)
		exp.xbeta.mat.prop	<- matrix(rep(exp.xbeta.prop, J), n, J)
		
		first.sum.prop 	<- colSums(exp.xbeta.mat.prop*ind.r_d)
		
		h.exp.xbeta.mat.prop 	<- - h.mat * exp.xbeta.mat.prop
		h.exp.xbeta.mat.prop[h.exp.xbeta.mat.prop > -10^(-7)] <- -10^(-7)
		second.sum.prop 	<- colSums(log(1 - exp(h.exp.xbeta.mat.prop))*ind.d)

		loglh.prop <- sum(-h*first.sum.prop + second.sum.prop)
		
		if(j == 1){
			logprior.prop   <- -1/2/sigmaSq*(be.prop[j]^2*(1/tauSq[j]+1/wSq[j])				-2*(be.prop[j]*be.ini[j+1]/wSq[j]))

			logprior.ini   <- -1/2/sigmaSq*(be.ini[j]^2*(1/tauSq[j]+1/wSq[j])				-2*(be.ini[j]*be.ini[j+1]/wSq[j]))
		}

		if(j == p){
			logprior.prop   <- -1/2/sigmaSq*(be.prop[j]^2*(1/tauSq[j]+1/wSq[j-1])					 -2*(be.ini[j-1]*be.prop[j]/wSq[j-1]))

			logprior.ini   <- -1/2/sigmaSq*(be.ini[j]^2*(1/tauSq[j]+1/wSq[j-1])					  -2*(be.ini[j-1]*be.ini[j]/wSq[j-1]))
		}


		if(j != 1 & j != p){
			logprior.prop   <- -1/2/sigmaSq*(be.prop[j]^2*(1/tauSq[j]+1/wSq[j-1]+1/								wSq[j])-2*(be.ini[j-1]*be.prop[j]/wSq[j-1]+								be.prop[j]*be.ini[j+1]/wSq[j]))

			logprior.ini   <- -1/2/sigmaSq*(be.ini[j]^2*(1/tauSq[j]+1/wSq[j-1]+1/wSq[j])					  -2*(be.ini[j-1]*be.ini[j]/wSq[j-1]+be.ini[j]*be.ini[j+1]/wSq[j]))
		}				
			

		logprop.prop    <- dnorm(be.prop[j] , mean = beta.prop.me[j], sd = sqrt(beta.prop.var), log = TRUE)
		logprop.ini     <- dnorm(be.ini[j]  , mean = beta.prop.me[j], sd = sqrt(beta.prop.var), log = TRUE)

		logR  <- loglh.prop - loglh.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop

		u = log(runif(1)) < logR
	
		if(u == 1){
			be.ini[j] <- be.prop[j]
			xbeta	<- xbeta.prop
			}

		accept[j]<-accept[j] + u
		
		} # end of for loop for j
		
	list(beta.ini = be.ini, accept = accept, xbeta = xbeta)
		
	}

