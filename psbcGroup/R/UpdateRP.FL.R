
UpdateRP.FL		<- function(survObj, priorPara, mcmcPara, ini){
	
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
		
		xbeta[xbeta > 700] <- 700
		exp.xbeta		<- exp(xbeta)
		x.exp.xbeta		<- x[,j] * exp.xbeta
		exp.xbeta.mat	<- matrix(rep(exp.xbeta, J), n, J)
		x.exp.xbeta.mat	<- matrix(rep(x.exp.xbeta, J), n, J)
		
		D1.1st <- -h * colSums(x.exp.xbeta.mat*ind.r_d)
		
		h.mat				<- matrix(rep(h, n), n, J, byrow = T)
		h.exp.xbeta.mat 	<- -h.mat * exp.xbeta.mat
		h.exp.xbeta.mat[h.exp.xbeta.mat > -10^(-7)] <- -10^(-7)
		exp.h.exp.xbeta.mat <- exp(h.exp.xbeta.mat)
		
		D1.2nd.den <- 1 - exp.h.exp.xbeta.mat
		D1.2nd.num <- exp.h.exp.xbeta.mat * x.exp.xbeta.mat
		D1.2nd <- h * colSums(D1.2nd.num/D1.2nd.den*ind.d)
		
		
		if(j == 1){
			D1.sd.be   <- -1/sigmaSq*(be.ini[j]*(1/tauSq[j]+1/wSq[j])-(be.ini[j+1]/wSq[j]))
		}

		if(j != 1 & j != p){
			D1.sd.be   <- -1/sigmaSq*(be.ini[j]*(1/tauSq[j]+1/wSq[j-1]+1/wSq[j])-(be.ini[j-1]/wSq[j-1]+be.ini[j+1]/wSq[j]))
		}		
		
		if(j == p){
			D1.sd.be   <- -1/sigmaSq*(be.ini[j]*(1/tauSq[j]+1/wSq[j-1])-(be.ini[j-1]/wSq[j-1]))
		}		
		
		
		D1 <- sum(D1.1st + D1.2nd) + D1.sd.be
		
		x.sq.exp.xbeta		<- x[,j]^2 * exp.xbeta
		x.sq.exp.xbeta.mat	<- matrix(rep(x.sq.exp.xbeta, J), n, J)
							
		D2.1st <- -h * colSums(x.sq.exp.xbeta.mat*ind.r_d)
				
		D2.2nd.den <- D1.2nd.den^2
		D2.2nd.num <- exp.h.exp.xbeta.mat * x.sq.exp.xbeta.mat * 							(1-exp.h.exp.xbeta.mat +h.exp.xbeta.mat)
		D2.2nd <- h * colSums(D2.2nd.num/D2.2nd.den*ind.d)
		
		if(j == 1){
			D2.sd.be   <- -1/sigmaSq*(1/tauSq[j]+1/wSq[j])
		}

		if(j != 1 & j != p){
			D2.sd.be   <- -1/sigmaSq*(1/tauSq[j]+1/wSq[j-1]+1/wSq[j])
		}		
		
		if(j == p){
			D2.sd.be   <- -1/sigmaSq*(1/tauSq[j]+1/wSq[j-1])
		}		
				
		
		D2 <- sum(D2.1st + D2.2nd) + D2.sd.be

		# proposal value
	
		be.prop.me  <- be.ini[j] - D1 / D2;
		be.prop.var <- -2.4^2 / D2;

		be.prop <- be.ini
		
		be.prop[j]	<- rnorm(1, mean = be.prop.me, sd = sqrt(be.prop.var));
				

		# Calculating acceptance probability

		xbeta.prop	<- xbeta - x[,j] * be.ini[j] + x[,j] * be.prop[j]
		
		xbeta.prop[xbeta.prop > 700] <- 700
		exp.xbeta.prop		<- exp(xbeta.prop)
		x.exp.xbeta.prop	<- x[,j] * exp.xbeta.prop
		exp.xbeta.prop.mat	<- matrix(rep(exp.xbeta.prop, J), n, J)
		x.exp.xbeta.prop.mat	<- matrix(rep(x.exp.xbeta.prop, J), n, J)
		
		D1.1st.prop <- -h * colSums(x.exp.xbeta.prop.mat*ind.r_d)
		
		h.exp.xbeta.prop.mat 	<- -h.mat * exp.xbeta.prop.mat
		h.exp.xbeta.prop.mat[h.exp.xbeta.prop.mat > -10^(-7)] <- -10^(-7)
		exp.h.exp.xbeta.prop.mat <- exp(h.exp.xbeta.prop.mat)
		
		D1.2nd.den.prop <- 1 - exp.h.exp.xbeta.prop.mat
		D1.2nd.num.prop <- exp.h.exp.xbeta.prop.mat * x.exp.xbeta.prop.mat
		D1.2nd.prop <- h * colSums(D1.2nd.num.prop/D1.2nd.den.prop*ind.d)

		if(j == 1){
			D1.sd.be.prop   <- -1/sigmaSq*(be.prop[j]*(1/tauSq[j]+1/wSq[j])-(be.ini[j+1]/wSq[j]))
		}

		if(j != 1 & j != p){
			D1.sd.be.prop   <- -1/sigmaSq*(be.prop[j]*(1/tauSq[j]+1/wSq[j-1]+1/wSq[j])-(be.ini[j-1]/wSq[j-1]+be.ini[j+1]/wSq[j]))
		}		
		
		if(j == p){
			D1.sd.be.prop   <- -1/sigmaSq*(be.prop[j]*(1/tauSq[j]+1/wSq[j-1])-(be.ini[j-1]/wSq[j-1]))
		}

		
		D1.prop <- sum(D1.1st.prop + D1.2nd.prop) + D1.sd.be.prop
		
		x.sq.exp.xbeta.prop		<- x[,j]^2 * exp.xbeta.prop
		x.sq.exp.xbeta.prop.mat	<- matrix(rep(x.sq.exp.xbeta.prop, J), n, J)
							
		D2.1st.prop <- -h * colSums(x.sq.exp.xbeta.prop.mat*ind.r_d)
				
		D2.2nd.den.prop <- D1.2nd.den.prop^2
		D2.2nd.num.prop <- exp.h.exp.xbeta.prop.mat * x.sq.exp.xbeta.prop.mat * 							(1-exp.h.exp.xbeta.prop.mat +h.exp.xbeta.prop.mat)
		D2.2nd.prop <- h * colSums(D2.2nd.num.prop/D2.2nd.den.prop*ind.d)
					
		D2.prop <- sum(D2.1st.prop + D2.2nd.prop) + D2.sd.be

		be.prop.me.ini  <- be.prop[j] - D1.prop / D2.prop;
		be.prop.var.ini <- -2.4^2 / D2.prop;

		#### -
				
		first.sum 	<- colSums(exp.xbeta.mat*ind.r_d)
		second.sum 	<- colSums(log(D1.2nd.den)*ind.d)

		loglh.ini <- sum(-h*first.sum + second.sum)
		
		first.sum.prop 	<- colSums(exp.xbeta.prop.mat*ind.r_d)
		second.sum.prop <- colSums(log(D1.2nd.den.prop)*ind.d)

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
			

		logprop.prop    <- dnorm(be.prop[j] , mean = be.prop.me.ini, sd = sqrt(be.prop.var.ini), log = TRUE)
		logprop.ini     <- dnorm(be.ini[j]  , mean = be.prop.me, sd = sqrt(be.prop.var), log = TRUE)

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













