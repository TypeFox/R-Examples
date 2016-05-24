psbcFL <-
function(survObj, priorPara, initial, rw = FALSE, mcmcPara, num.reps, thin, chain = 1, save = 1000){
	
	survObj$n <- n <- length(survObj$t)
	survObj$p <- p <- dim(survObj$x)[2]	
	
	eta0	<- priorPara$eta0
	kappa0	<- priorPara$kappa0
	c0		<- priorPara$c0
	r1		<- priorPara$r1
	r2		<- priorPara$r2
	delta1	<- priorPara$delta1
	delta2	<- priorPara$delta2	
	s		<- priorPara$s
	J		<- priorPara$J	<- length(priorPara$s)
	
	intv 				<- setting.interval(survObj$t, survObj$di, 							priorPara$s, priorPara$J)
	priorPara$ind.r		<- intv$ind.r
	priorPara$ind.d		<- intv$ind.d
	priorPara$ind.r_d	<- intv$ind.r_d	
	priorPara$d			<- intv$d
	
	ini	<- initial
	
	beta.ini			<- ini$beta.ini
	lambda1Sq			<- ini$lambda1Sq
	lambda2Sq			<- ini$lambda2Sq
	sigmaSq				<- ini$sigmaSq
	tauSq				<- ini$tauSq
	h					<- ini$h	
	wSq					<- ini$wSq

	mcmcPara$beta.prop.me	<- beta.ini	
	
	covInv <- matrix(0, p, p)
	diag(covInv) <- 1/tauSq + c(0, 1/wSq) + c(1/wSq, 0)
	covInv[1,2] <- -1/wSq[1]
	covInv[p,(p-1)] <- -1/wSq[p-1]

	for(j in 2:(p-1)){
		covInv[j,j-1] <- covInv[j-1,j] <- -1/wSq[j-1]
	}	
	
	ini$covInv	<- covInv	
	ini$xbeta	<- as.vector(survObj$x %*% beta.ini)
	
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
	wSq.p						<- wSq	
	mcmcOutcome$sigmaSq.p		<- sigmaSq
	mcmcOutcome$lambda1Sq.p		<- lambda1Sq
	mcmcOutcome$lambda2Sq.p		<- lambda2Sq
		
	mcmcOutcome$accept.beta		<- c(rep(0, p))
	
	
	outcomeSum	<- list()
	
	dir.create('mcmcOutcome', showWarnings = FALSE)
		

	# MCMC sampling

	for(M in 1:num.reps){	
	
		cat("Chain", chain, "Iteration", M, fill=TRUE);

		
		# Updating regression parameters
		
		if(rw == FALSE){
			sampleRP	<- UpdateRP.FL(survObj, priorPara, mcmcPara, ini)
			}
		if(rw == TRUE){
			sampleRP	<- UpdateRP.FLrw(survObj, priorPara, mcmcPara, ini)
			}		
		
		beta.ini	<- ini$beta.ini	<- sampleRP$beta.ini
		xbeta		<- ini$xbeta	<- sampleRP$xbeta
		mcmcOutcome$accept.beta	<- mcmcOutcome$accept.beta + sampleRP$accept
		
		
		# Updating increments in cumulative hazards
				
		h	<- ini$h	<- UpdateBH(survObj, priorPara, ini)
	
	
		# Updating 1/tauSq
		
		tauSq	<- ini$tauSq	<- UpdateTau(survObj, priorPara, ini)
		
		
		# Updating 1/wSq
		
		wSq	<- ini$wSq	<- UpdateW.FL(survObj, priorPara, ini)
		
		
		# Updating Sigma_beta based on new tauSq and wSq

		covInv <- matrix(0, p, p)
		diag(covInv) <- 1/tauSq + c(0, 1/wSq) + c(1/wSq, 0)
		covInv[1,2] <- -1/wSq[1]
		covInv[p,(p-1)] <- -1/wSq[p-1]

		for(j in 2:(p-1)){
			covInv[j,j-1] <- covInv[j-1,j] <- -1/wSq[j-1]
		}	
	
		ini$covInv	<- covInv	
		
						
		# Updating sigmaSq
		
		sigmaSq	<- ini$sigmaSq	<- UpdateSigma.FL(survObj, priorPara, ini)
		

		# Updating lambda1Sq		
		
		lambda1Sq	<- ini$lambda1Sq	<- UpdateLambda1.FL(survObj, priorPara, ini)
		
		
		# Updating lambda2Sq		
		
		lambda2Sq	<- ini$lambda2Sq	<- UpdateLambda2.FL(survObj, priorPara, ini)
				

		###### storing posterior samples
	
	if(M %% thin == 0){	 		
		beta.p						<- rbind(beta.p, beta.ini, deparse.level = 0)
		h.p							<- rbind(h.p, h, deparse.level = 0)
		tauSq.p						<- rbind(tauSq.p, tauSq, deparse.level = 0)
		wSq.p						<- rbind(wSq.p, wSq, deparse.level = 0)
		mcmcOutcome$sigmaSq.p		<- c(mcmcOutcome$sigmaSq.p, sigmaSq)
		mcmcOutcome$lambda1Sq.p		<- c(mcmcOutcome$lambda1Sq.p, lambda1Sq)
		mcmcOutcome$lambda2Sq.p		<- c(mcmcOutcome$lambda2Sq.p, lambda2Sq)
		
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
			save(wSq.p, file = paste("mcmcOutcome/wSqAll.ch", chain, ".Rdata", sep = ""))
			save(h.p, file = paste("mcmcOutcome/hAll.ch", chain, ".Rdata", sep = ""))
			}
				
		} # the end of MCMC sampling
		
	
	list(beta.p = beta.p, h.p = h.p, tauSq.p = tauSq.p, mcmcOutcome = mcmcOutcome)
	
	} # end of "psbcGrp" function

