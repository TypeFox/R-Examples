psbcGL <-
function(survObj, priorPara, initial, rw = FALSE, mcmcPara, num.reps, thin, chain = 1, save = 1000){
	
	survObj$n <- n <- length(survObj$t)
	survObj$p <- p <- dim(survObj$x)[2]	
	
	eta0	<- priorPara$eta0
	kappa0	<- priorPara$kappa0
	c0		<- priorPara$c0
	r		<- priorPara$r
	delta	<- priorPara$delta	
	s		<- priorPara$s
	J		<- priorPara$J	<- length(priorPara$s)
	groupInd 	<- priorPara$groupInd
	groupNo		<- priorPara$groupNo <- unique(priorPara$groupInd)
	K			<- priorPara$K <- length(groupNo)
	m_k			<- priorPara$m_k
	
	m_k 	<- rep(NA, K)		
	for(i in 1:K){
		m_k[i] <- sum(groupInd == groupNo[i])
	}	
	priorPara$m_k	<- m_k
	
	intv 				<- setting.interval(survObj$t, survObj$di, 							priorPara$s, priorPara$J)
	priorPara$ind.r		<- intv$ind.r
	priorPara$ind.d		<- intv$ind.d
	priorPara$ind.r_d	<- intv$ind.r_d	
	priorPara$d			<- intv$d
	
	ini	<- initial
	
	beta.ini			<- ini$beta.ini
	lambdaSq			<- ini$lambdaSq
	sigmaSq				<- ini$sigmaSq
	tauSq				<- ini$tauSq
	h					<- ini$h	
	
	mcmcPara$beta.prop.me	<- beta.ini	

	tauSq.exp <- rep(NA, p)

	for(i in 1:K){
	tauSq.exp[groupInd == groupNo[i]] <- tauSq[i]
	}

	ini$sd.be      <- sqrt(sigmaSq*tauSq.exp)	
	ini$xbeta	<- as.vector(survObj$x %*% beta.ini)
	
	be.normSq <- c()

	for(i in 1:K){
		be.normSq[i] <- sum(beta.ini[which(groupInd == groupNo[i])]^2)
		}
		
	ini$be.normSq	<- be.normSq			
		
	H.star   <- alpha0  <- c()

	for (j in 1:J){
		H.star[j]	<- eta0 * s[j]^kappa0
		alpha0[j]	<- c0 * H.star[j] 
		}

	priorPara$hPriorSh	<- diff(c(0, alpha0))
	
			
	## for posterior samples
	
	mcmcOutcome	<- list()
	
	mcmcOutcome$initial			<- initial
	mcmcOutcome$priorPara		<- priorPara
	
	beta.p						<- beta.ini
	h.p							<- h
	tauSq.p						<- tauSq
	mcmcOutcome$sigmaSq.p		<- sigmaSq
	mcmcOutcome$lambdaSq.p		<- lambdaSq		
	mcmcOutcome$accept.beta		<- c(rep(0, p))
	
	
	outcomeSum	<- list()
	
	dir.create('mcmcOutcome', showWarnings = FALSE)
		

	# MCMC sampling

	for(M in 1:num.reps){	
	
		cat("Chain", chain, "Iteration", M, fill=TRUE);
		
		# Updating regression parameters
		
		if(rw == FALSE){
			sampleRP	<- UpdateRP(survObj, priorPara, mcmcPara, ini)
			}
		if(rw == TRUE){
			sampleRP	<- UpdateRPrw(survObj, priorPara, mcmcPara, ini)
			}
		beta.ini	<- ini$beta.ini	<- sampleRP$beta.ini
		xbeta		<- ini$xbeta	<- sampleRP$xbeta
		mcmcOutcome$accept.beta	<- mcmcOutcome$accept.beta + sampleRP$accept
		
		# Updating the squared norm of beta using new beta values
		
		for(i in 1:K){
		be.normSq <- ini$be.normSq[i] <- sum(beta.ini[which(groupInd == groupNo[i])]^2)
		}		
		
		# Updating increments in cumulative hazards
				
		h	<- ini$h	<- UpdateBH(survObj, priorPara, ini)
	
	
		# Updating 1/tauSq
		
		tauSq	<- ini$tauSq	<- UpdateTau.GL(survObj, priorPara, ini)
		
		# Updating tauSq.exp with new tauSq

		for(i in 1:K){
			tauSq.exp[groupInd == groupNo[i]] <- tauSq[i]
			}
							
		# Updating sigmaSq
		
		sigmaSq	<- ini$sigmaSq	<- UpdateSigma.GL(survObj, priorPara, ini)
		

		# Updating lambdaSq		
		
		lambdaSq	<- ini$lambdaSq	<- UpdateLambda.GL(survObj, priorPara, ini)
		
		##########################
		
		ini$sd.be      <- sqrt(sigmaSq*tauSq.exp)
			
		###### storing posterior samples
	
	if(M %% thin == 0){	 		
		beta.p						<- rbind(beta.p, beta.ini, deparse.level = 0)
		h.p							<- rbind(h.p, h, deparse.level = 0)
		tauSq.p						<- rbind(tauSq.p, tauSq, deparse.level = 0)
		mcmcOutcome$sigmaSq.p		<- c(mcmcOutcome$sigmaSq.p, sigmaSq)
		mcmcOutcome$lambdaSq.p		<- c(mcmcOutcome$lambdaSq.p, lambdaSq)
		
		mcmcOutcome$ini				<- ini
		
		}
			
		###### Tuning algorithm for the mean of the proposal density ###
		

		for(j in 1:survObj$p){
			if(M%/%thin > (20%/%thin)){
				if(beta.ini[j] == beta.p[(M%/%thin + 1 -(20%/%thin)),j]){    
					mcmcPara$beta.prop.me[j] <- beta.p[(M%/%thin + 1),j];
					}
				}
			}

		
		# saving the mcmc outcomes

		if(M %% save == 0 | M == num.reps){
			save(mcmcOutcome, file = paste("mcmcOutcome/otherAll.ch", chain, ".Rdata", sep = ""))
			save(beta.p, file = paste("mcmcOutcome/betaAll.ch", chain, ".Rdata", sep = ""))
			save(tauSq.p, file = paste("mcmcOutcome/tauSqAll.ch", chain, ".Rdata", sep = ""))
			save(h.p, file = paste("mcmcOutcome/hAll.ch", chain, ".Rdata", sep = ""))
			}
				
		} # the end of MCMC sampling
		
	
	list(beta.p = beta.p, h.p = h.p, tauSq.p = tauSq.p, mcmcOutcome = mcmcOutcome)
	
	} # end of "psbcGrp" function

