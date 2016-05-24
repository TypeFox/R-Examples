
######################################################################################################################################
######################################################################################################################################
### HiSSE -- BiSSE with hidden states
######################################################################################################################################
######################################################################################################################################

hisse <- function(phy, data, f=c(1,1), hidden.states=TRUE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=NULL, turnover.beta=c(0,0,0,0), eps.beta=c(0,0,0,0), timeslice=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, output.type="turnover", sann=FALSE, sann.its=10000, max.tol=.Machine$double.eps^.25){
	if(!is.null(root.p)) {
		root.type="user"
		root.p <- root.p / sum(root.p)	
		if(hidden.states ==TRUE & length(root.p)==2){
			root.p <- rep(root.p, 2)
			root.p <- root.p / sum(root.p)	
			warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
		}
	}
	
	if(is.null(trans.rate)){
		stop("Rate matrix needed. See TransMatMaker() to create one.")
	}
	
	if(hidden.states == TRUE & dim(trans.rate)[1]<4){
		stop("You chose a hidden state but this is not reflected in the transition matrix")
	}
	
	if(sum(turnover.beta) > 0 | sum(eps.beta) > 0){
		warning("You chose a time dependent model. This is currently untested -- Good luck.")
	}

	if(!is.null(timeslice)){
		stop("You chose a time slice model. This is currently unavailable.")
	}
	
	#Some basic checks:
	if(!length(turnover.anc) == 4){
		turnover.anc<-c(turnover.anc,0,0)
	}
	if(!length(eps.anc) == 4){
		eps.anc<-c(eps.anc,0,0)
	}
	if(!length(turnover.beta) == 4){
		eps.anc<-c(eps.anc,0,0)
	}
	if(!length(eps.beta) == 4){
		eps.anc<-c(eps.anc,0,0)
	}
	phy$node.label <- NULL
	#Sets the main parameters to be used in the model:
	if(sum(eps.anc)==0){
		eps.anc=c(0,0,0,0)
		eps.beta=c(0,0,0,0)
	}
	
	#Initialize our pars vector and build from there:
	pars = c(turnover.anc)
	#Add in extinction fraction:
	eps.anc.tmp = eps.anc
	eps.anc.tmp[which(eps.anc.tmp>0)] = (eps.anc.tmp[which(eps.anc.tmp>0)] + max(pars))
	pars = c(pars, eps.anc.tmp)
	#Add in transition rates:
	if(hidden.states == FALSE){
		trans.rate.tmp = c(trans.rate[!is.na(trans.rate)][1], 0, 0, trans.rate[!is.na(trans.rate)][2], rep(0, 8))
	}else{
		trans.rate.tmp = trans.rate[!is.na(trans.rate)]
	}
	trans.rate.tmp[which(trans.rate.tmp > 0)] = (trans.rate.tmp[which(trans.rate.tmp > 0)] + max(pars))
	pars = c(pars, trans.rate.tmp)
	#Add in turnover trend parameters
	turnover.beta.alpha = turnover.beta
	turnover.beta.beta = turnover.beta
 	turnover.beta.beta[turnover.beta.beta>0] = turnover.beta.beta[turnover.beta.beta > 0] + max(turnover.beta.beta)
	turnover.beta.tmp <- c(turnover.beta.alpha, turnover.beta.beta)
	turnover.beta.tmp[which(turnover.beta.tmp > 0)] = (turnover.beta.tmp[which(turnover.beta.tmp > 0)] + max(pars))
	pars = c(pars, turnover.beta.tmp)		
	#Add in eps trend parameters:
	eps.beta.alpha = eps.beta
	eps.beta.beta = eps.beta
 	eps.beta.beta[eps.beta.beta>0] = eps.beta.beta[eps.beta.beta > 0] + max(eps.beta.beta)
	eps.beta.tmp <- c(eps.beta.alpha, eps.beta.beta)
	eps.beta.tmp[which(eps.beta.tmp > 0)] = (eps.beta.tmp[which(eps.beta.tmp > 0)] + max(pars))
	pars = c(pars, eps.beta.tmp)
	if(is.null(timeslice)){
		turnover.slice = c(0,0,0,0)
		eps.slice = c(0,0,0,0)
		trans.rate.slice = rep(0, 12)	
	}else{
		turnover.slice = turnover.anc
		eps.slice = eps.anc
		trans.rate.slice = trans.rate	
	}	
	#Add in turnover slice factors:
	turnover.slice.tmp = turnover.slice
	turnover.slice.tmp[which(turnover.slice.tmp > 0)] = (turnover.slice.tmp[which(turnover.slice.tmp > 0)] + max(pars))
	pars = c(pars, turnover.slice.tmp)		
	#Add in eps slice factors:
	eps.slice.tmp = eps.slice
	eps.slice.tmp[which(eps.slice.tmp > 0)] = (eps.slice.tmp[which(eps.slice.tmp > 0)] + max(pars))
	pars = c(pars, eps.slice.tmp)		
	#Add in transition rate slice factors:
	trans.rate.slice.tmp = trans.rate.slice
	trans.rate.slice.tmp[which(trans.rate.slice.tmp > 0)] = (trans.rate.slice.tmp[which(trans.rate.slice.tmp > 0)] + max(pars))
	pars = c(pars, trans.rate.slice.tmp)		
	
	np <- max(pars)
	pars[pars==0] <- np+1
	
	cat("Initializing...", "\n")
	
	data.new<-data.frame(data[,2], data[,2], row.names=data[,1])
	data.new<-data.new[phy$tip.label,]
	#This is used to scale starting values to account for sampling:
	if(length(f) == 2){
		samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f)
	}else{
		if(length(f) == Ntip(phy)){
			samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / mean(f))
		}else{
			stop("The vector of sampling frequencies does not match the number of tips in the tree.")
		}
	}

	if(sum(eps.anc)==0){
		init.pars <- starting.point.generator(phy, 2, samp.freq.tree, yule=TRUE)
		names(init.pars) <- NULL
		def.set.pars <- c(rep(log(init.pars[1]+init.pars[3]), 4), rep(log(init.pars[3]/init.pars[1]),4), rep(log(init.pars[5]), 12), rep(log(1), 36))
	}else{
		init.pars <- starting.point.generator(phy, 2, samp.freq.tree, yule=FALSE)
		names(init.pars) <- NULL
		init.eps = init.pars[3]/init.pars[1]
		if(init.eps == 0){
			init.eps = 1e-6
		}
		def.set.pars <- c(rep(log(init.pars[1]+init.pars[3]), 4), rep(log(init.eps),4), rep(log(init.pars[5]), 12), rep(log(1), 36))
	}
	#Set initials using estimates from constant bd model:
	np.sequence <- 1:np
	ip <- numeric(np)
	for(i in np.sequence){
		ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
	}

	lower <- rep(-20, length(ip))
	upper <- -lower
	
	if(sann == FALSE){
		cat("Finished. Beginning subplex routine...", "\n")
		out = subplex(ip, fn=DevOptimize, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, timeslice=timeslice, np=np)
	}else{
		cat("Finished. Beginning simulated annealing...", "\n")
		out.sann = GenSA(ip, fn=DevOptimize, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, timeslice=timeslice, np=np)
		cat("Finished. Refining using subplex routine...", "\n")
		out = subplex(out.sann$par, fn=DevOptimize, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, timeslice=timeslice, np=np)
	}
	solution <- numeric(length(pars))
#	solution[] <- c(exp(out$solution), 0)[pars]
	solution[] <- c(exp(out$par), 0)[pars]
#	loglik = -out$objective
	loglik = -out$value
	
	cat("Finished. Summarizing results...", "\n")

	#Some cleanup to make the output look pretty:
	solution.tmp = solution[21:56]
	solution.tmp[solution.tmp==0] = 1
	solution[21:56] = solution.tmp
	
	obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, timeslice=timeslice, phy=phy, data=data, output.type=output.type, max.tol=max.tol) 
	class(obj) = "hisse.fit"		
	
	return(obj)		
}



######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################

#Function used for optimizing parameters:
DevOptimize <- function(p, pars, phy, data, f, hidden.states, condition.on.survival, root.type, root.p, timeslice, np) {
	#Generates the final vector with the appropriate parameter estimates in the right place:
	p.new <- exp(p)
	model.vec <- numeric(length(pars))
	model.vec[] <- c(p.new, 0)[pars]
	model.vec.tmp = model.vec[21:56]
	model.vec.tmp[model.vec.tmp==0] = 1
	model.vec[21:56] = model.vec.tmp
	
	cache = ParametersToPass(phy, data, model.vec, f=f, timeslice=timeslice, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
	cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
	cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])

	cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
	cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
	cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
	cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])
	logl <- DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
	return(-logl)
}


#Taken from the BiSSE code -- credit goes to Rich FitzJohn:
starting.point.tree <- function(phy, yule=FALSE) {
	p.yule <- c(yule(phy)$lambda, 0)
	if(yule){
		p.yule
	}else{
		suppressWarnings(c(birthdeath(phy)$para[2] / (1-birthdeath(phy)$para[1]), ((birthdeath(phy)$para[1] * birthdeath(phy)$para[2]) / (1-birthdeath(phy)$para[1]))))
	}
}


starting.point.generator <- function(phy, k, samp.freq.tree, q.div=5, yule=FALSE) {
	pars.bd <- suppressWarnings(starting.point.tree(phy, yule))
	#Rescale parameters to account for sampling, if necessary, using Stadler 2013: 
	pars.bd[1] = pars.bd[1] / samp.freq.tree
	pars.bd[2] = pars.bd[2] - (pars.bd[1]*samp.freq.tree) * (1 - 1/samp.freq.tree)
	r <- if  ( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]
	p <- rep(c(pars.bd, r / q.div), c(k, k, k * (k-1)))
	names(p) <- NULL
	p
}


######################################################################################################################################
######################################################################################################################################
### The down pass that carries out the integration and returns the likelihood: 
######################################################################################################################################
######################################################################################################################################

DownPass <- function(phy, cache, hidden.states, bad.likelihood=-10000000000, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL) {
	#Some preliminaries:
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	phy <- reorder(phy, "pruningwise")
	anc <- unique(phy$edge[,1])		
	TIPS <- 1:nb.tip

	if(hidden.states == FALSE){
		compD <- matrix(0, nrow=nb.tip + nb.node, ncol=2)
		compE <- matrix(0, nrow=nb.tip + nb.node, ncol=2)
	}else{
		compD <- matrix(0, nrow=nb.tip + nb.node, ncol=4)
		compE <- matrix(0, nrow=nb.tip + nb.node, ncol=4)		
	}
		
	#Initializes the tip sampling and sets internal nodes to be zero:
	ncols = dim(compD)[2]
	if(length(cache$f) == 2){
		for(i in 1:(nb.tip)){
			compD[i,] <- cache$f * cache$states[i,]
			compE[i,] <- rep((1-cache$f), ncols/2)
		}
	}else{
		for(i in 1:(nb.tip)){
			compD[i,] <- cache$f[i] * cache$states[i,]
			compE[i,] <- rep((1-cache$f[i]), ncols/2)
		}
	}
	logcomp <- c()
	#Start the postorder traversal indexing lists by node number: 
	for (i in seq(from = 1, length.out = nb.node)) {
		#A vector of all the internal nodes:
		focal <- anc[i]
		desRows <- which(phy$edge[,1]==focal)
		desNodes <- phy$edge[desRows,2]
		#Note: when the tree has been reordered branching.times are no longer valid. Fortunately, we extract this information in the initial cache setup. Also focal is the rootward node, whereas desNodes represent a vector of all descendant nodes:
		cache$rootward.age <- cache$split.times[which(names(cache$split.times)==focal)]
		
		if(hidden.states == FALSE){
			v = c(1,1)
		}else{
			v = c(1,1,1,1)
		}
		phi <- c()
		for (desIndex in sequence(length(desRows))){
			cache$focal.edge.length <- phy$edge.length[desRows[desIndex]]
			cache$tipward.age <- cache$rootward.age - cache$focal.edge.length
			#Strange rounding errors. A tip age should be zero. This ensures that:
			if(cache$tipward.age < .Machine$double.eps^0.5){
				cache$tipward.age = 0
			}
			cache$node.D <- compD[desNodes[desIndex],]
			cache$node.E <- compE[desNodes[desIndex],]
			##Call to lsoda that utilizes C code. Requires a lot of inputs. Note that for now we hardcode the NUMELEMENTS arguments. The reason for this is because with lsoda we can only pass a vector of parameters.
			if(hidden.states == FALSE){
				pars <- list(cache$tot_time, cache$timeslice, cache$turnover.trend.alpha0, cache$turnover.trend.beta0, cache$turnover.beta.factor0, cache$turnover.slice.factor0, cache$eps.trend.alpha0, cache$eps.trend.beta0, cache$eps.beta.factor0, cache$eps.slice.factor0, cache$turnover.trend.alpha1, cache$turnover.trend.beta1, cache$turnover.beta.factor1, cache$turnover.slice.factor1, cache$eps.trend.alpha1, cache$eps.trend.beta1, cache$eps.beta.factor1, cache$eps.slice.factor1, cache$x_turnover0, cache$x_eps0, cache$x_turnover1, cache$x_eps1, cache$q01, cache$q10, cache$q01_slice.factor, cache$q10_slice.factor, cache$focal.edge.length, cache$tipward.age)
				NUMELEMENTS <- 28 #needed for passing in vector to C
				padded.pars <- rep(0, NUMELEMENTS)
				pars <- c(unlist(pars))
				stopifnot(length(padded.pars)<=NUMELEMENTS)
				padded.pars[sequence(length(pars))]<-pars
				yini <-c(E0=cache$node.E[1], E1=cache$node.E[2], D0=cache$node.D[1], D1=cache$node.D[2])
				times=c(cache$tipward.age, cache$rootward.age)							
				prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_bisse", padded.pars, initfunc="initmod_bisse", dllname = "hisse", rtol=1e-8, atol=1e-8, hmax=10)
			}else{
				pars <- list(cache$tot_time, cache$timeslice, cache$turnover.trend.alpha0, cache$turnover.trend.beta0, cache$turnover.beta.factor0, cache$turnover.slice.factor0, cache$eps.trend.alpha0, cache$eps.trend.beta0, cache$eps.beta.factor0, cache$eps.slice.factor0, cache$turnover.trend.alpha1, cache$turnover.trend.beta1, cache$turnover.beta.factor1, cache$turnover.slice.factor1, cache$eps.trend.alpha1, cache$eps.trend.beta1, cache$eps.beta.factor1, cache$eps.slice.factor1, cache$turnover.trend.alphaA, cache$turnover.trend.betaA, cache$turnover.beta.factorA, cache$turnover.slice.factorA, cache$eps.trend.alphaA, cache$eps.trend.betaA, cache$eps.beta.factorA, cache$eps.slice.factorA, cache$turnover.trend.alphaB, cache$turnover.trend.betaB, cache$turnover.beta.factorB, cache$turnover.slice.factorB, cache$eps.trend.alphaB, cache$eps.trend.betaB, cache$eps.beta.factorB, cache$eps.slice.factorB, cache$x_turnover0, cache$x_eps0, cache$x_turnover1, cache$x_eps1, cache$x_turnoverA, cache$x_epsA, cache$x_turnoverB, cache$x_epsB, cache$q01, cache$q10, cache$q0A, cache$qA0, cache$q1B, cache$qB1, cache$q0B, cache$qB0, cache$q1A, cache$qA1, cache$qBA, cache$qAB, cache$q01_slice.factor, cache$q10_slice.factor, cache$q0A_slice.factor, cache$qA0_slice.factor, cache$q1B_slice.factor, cache$qB1_slice.factor, cache$q0B_slice.factor, cache$qB0_slice.factor, cache$q1A_slice.factor, cache$qA1_slice.factor, cache$qBA_slice.factor, cache$qAB_slice.factor, cache$focal_edge_length, cache$tipward_age)
				NUMELEMENTS <- 68 #needed for passing in vector to C
				padded.pars <- rep(0, NUMELEMENTS)
				pars <- c(unlist(pars))
				stopifnot(length(padded.pars)<=NUMELEMENTS)
				padded.pars[sequence(length(pars))]<-pars				
				yini <-c(E0=cache$node.E[1], E1=cache$node.E[2], EA=cache$node.E[3], EB=cache$node.E[4], D0=cache$node.D[1], D1=cache$node.D[2], DA=cache$node.D[3], DB=cache$node.D[4])
				times=c(cache$tipward.age, cache$rootward.age)							
				prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_hisse", padded.pars, initfunc="initmod_hisse", dllname = "hisse", rtol=1e-8, atol=1e-8, hmax=10)
			}

			######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
			if(attributes(prob.subtree.cal.full)$istate[1] < 0){
				return(bad.likelihood)
			}else{
				prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
			}
			##############################################################################
			if(hidden.states == FALSE){
				if(is.nan(prob.subtree.cal[3]) | is.nan(prob.subtree.cal[4])){
					return(bad.likelihood)
				}										
				if(prob.subtree.cal[3]<0 | prob.subtree.cal[4]<0){
					return(bad.likelihood)
				}
				#Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
				phi <- c(phi, prob.subtree.cal[1:2])
				v <- v * prob.subtree.cal[3:4]
			}else{
				if(is.nan(prob.subtree.cal[5]) | is.nan(prob.subtree.cal[6]) | is.nan(prob.subtree.cal[7]) | is.nan(prob.subtree.cal[8])){
					return(bad.likelihood)
				}										
				if(prob.subtree.cal[5]<0 | prob.subtree.cal[6]<0 | prob.subtree.cal[7]<0 | prob.subtree.cal[8]<0){
					return(bad.likelihood)
				}
				#Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
				phi <- c(phi, prob.subtree.cal[1:4])
				v <- v * prob.subtree.cal[5:8]				
			}
		}
		#C call to set_birth_void -- NOTE: The first input is zero as we need to declare the birth_rate. It gets written over and is now the first element in the list that is returned. Everything else should be self explanatory.
		if(hidden.states == TRUE){			
			lambda0 <- .C("set_birth_hisse_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$timeslice), as.double(cache$turnover.trend.alpha0), as.double(cache$turnover.trend.beta0), as.double(cache$turnover.beta.factor0), as.double(cache$turnover.slice.factor0), as.double(cache$eps.trend.alpha0), as.double(cache$eps.trend.beta0), as.double(cache$eps.beta.factor0), as.double(cache$eps.slice.factor0), as.double(cache$turnover.trend.alpha1), as.double(cache$turnover.trend.beta1), as.double(cache$turnover.beta.factor1), as.double(cache$turnover.slice.factor1), as.double(cache$eps.trend.alpha1), as.double(cache$eps.trend.beta1), as.double(cache$eps.beta.factor1), as.double(cache$eps.slice.factor1), as.double(cache$turnover.trend.alphaA), as.double(cache$turnover.trend.betaA), as.double(cache$turnover.beta.factorA), as.double(cache$turnover.slice.factorA), as.double(cache$eps.trend.alphaA), as.double(cache$eps.trend.betaA), as.double(cache$eps.beta.factorA), as.double(cache$eps.slice.factorA), as.double(cache$turnover.trend.alphaB), as.double(cache$turnover.trend.betaB), as.double(cache$turnover.beta.factorB), as.double(cache$turnover.slice.factorB), as.double(cache$eps.trend.alphaB), as.double(cache$eps.trend.betaB), as.double(cache$eps.beta.factorB), as.double(cache$eps.slice.factorB), as.double(cache$x_turnover0), as.double(cache$x_eps0), as.double(cache$x_turnover1), as.double(cache$x_eps1), as.double(cache$x_turnoverA), as.double(cache$x_epsA), as.double(cache$x_turnoverB), as.double(cache$x_epsB), as.double(cache$q01), as.double(cache$q10), as.double(cache$q0A), as.double(cache$qA0), as.double(cache$q1B), as.double(cache$qB1), as.double(cache$q0B), as.double(cache$qB0), as.double(cache$q1A), as.double(cache$qA1), as.double(cache$qBA), as.double(cache$qAB), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(0))
			lambda1 <- .C("set_birth_hisse_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$timeslice), as.double(cache$turnover.trend.alpha0), as.double(cache$turnover.trend.beta0), as.double(cache$turnover.beta.factor0), as.double(cache$turnover.slice.factor0), as.double(cache$eps.trend.alpha0), as.double(cache$eps.trend.beta0), as.double(cache$eps.beta.factor0), as.double(cache$eps.slice.factor0), as.double(cache$turnover.trend.alpha1), as.double(cache$turnover.trend.beta1), as.double(cache$turnover.beta.factor1), as.double(cache$turnover.slice.factor1), as.double(cache$eps.trend.alpha1), as.double(cache$eps.trend.beta1), as.double(cache$eps.beta.factor1), as.double(cache$eps.slice.factor1), as.double(cache$turnover.trend.alphaA), as.double(cache$turnover.trend.betaA), as.double(cache$turnover.beta.factorA), as.double(cache$turnover.slice.factorA), as.double(cache$eps.trend.alphaA), as.double(cache$eps.trend.betaA), as.double(cache$eps.beta.factorA), as.double(cache$eps.slice.factorA), as.double(cache$turnover.trend.alphaB), as.double(cache$turnover.trend.betaB), as.double(cache$turnover.beta.factorB), as.double(cache$turnover.slice.factorB), as.double(cache$eps.trend.alphaB), as.double(cache$eps.trend.betaB), as.double(cache$eps.beta.factorB), as.double(cache$eps.slice.factorB), as.double(cache$x_turnover0), as.double(cache$x_eps0), as.double(cache$x_turnover1), as.double(cache$x_eps1), as.double(cache$x_turnoverA), as.double(cache$x_epsA), as.double(cache$x_turnoverB), as.double(cache$x_epsB), as.double(cache$q01), as.double(cache$q10), as.double(cache$q0A), as.double(cache$qA0), as.double(cache$q1B), as.double(cache$qB1), as.double(cache$q0B), as.double(cache$qB0), as.double(cache$q1A), as.double(cache$qA1), as.double(cache$qBA), as.double(cache$qAB), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(1))
			lambdaA <- .C("set_birth_hisse_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$timeslice), as.double(cache$turnover.trend.alpha0), as.double(cache$turnover.trend.beta0), as.double(cache$turnover.beta.factor0), as.double(cache$turnover.slice.factor0), as.double(cache$eps.trend.alpha0), as.double(cache$eps.trend.beta0), as.double(cache$eps.beta.factor0), as.double(cache$eps.slice.factor0), as.double(cache$turnover.trend.alpha1), as.double(cache$turnover.trend.beta1), as.double(cache$turnover.beta.factor1), as.double(cache$turnover.slice.factor1), as.double(cache$eps.trend.alpha1), as.double(cache$eps.trend.beta1), as.double(cache$eps.beta.factor1), as.double(cache$eps.slice.factor1), as.double(cache$turnover.trend.alphaA), as.double(cache$turnover.trend.betaA), as.double(cache$turnover.beta.factorA), as.double(cache$turnover.slice.factorA), as.double(cache$eps.trend.alphaA), as.double(cache$eps.trend.betaA), as.double(cache$eps.beta.factorA), as.double(cache$eps.slice.factorA), as.double(cache$turnover.trend.alphaB), as.double(cache$turnover.trend.betaB), as.double(cache$turnover.beta.factorB), as.double(cache$turnover.slice.factorB), as.double(cache$eps.trend.alphaB), as.double(cache$eps.trend.betaB), as.double(cache$eps.beta.factorB), as.double(cache$eps.slice.factorB), as.double(cache$x_turnover0), as.double(cache$x_eps0), as.double(cache$x_turnover1), as.double(cache$x_eps1), as.double(cache$x_turnoverA), as.double(cache$x_epsA), as.double(cache$x_turnoverB), as.double(cache$x_epsB), as.double(cache$q01), as.double(cache$q10), as.double(cache$q0A), as.double(cache$qA0), as.double(cache$q1B), as.double(cache$qB1), as.double(cache$q0B), as.double(cache$qB0), as.double(cache$q1A), as.double(cache$qA1), as.double(cache$qBA), as.double(cache$qAB), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(2))
			lambdaB <- .C("set_birth_hisse_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$timeslice), as.double(cache$turnover.trend.alpha0), as.double(cache$turnover.trend.beta0), as.double(cache$turnover.beta.factor0), as.double(cache$turnover.slice.factor0), as.double(cache$eps.trend.alpha0), as.double(cache$eps.trend.beta0), as.double(cache$eps.beta.factor0), as.double(cache$eps.slice.factor0), as.double(cache$turnover.trend.alpha1), as.double(cache$turnover.trend.beta1), as.double(cache$turnover.beta.factor1), as.double(cache$turnover.slice.factor1), as.double(cache$eps.trend.alpha1), as.double(cache$eps.trend.beta1), as.double(cache$eps.beta.factor1), as.double(cache$eps.slice.factor1), as.double(cache$turnover.trend.alphaA), as.double(cache$turnover.trend.betaA), as.double(cache$turnover.beta.factorA), as.double(cache$turnover.slice.factorA), as.double(cache$eps.trend.alphaA), as.double(cache$eps.trend.betaA), as.double(cache$eps.beta.factorA), as.double(cache$eps.slice.factorA), as.double(cache$turnover.trend.alphaB), as.double(cache$turnover.trend.betaB), as.double(cache$turnover.beta.factorB), as.double(cache$turnover.slice.factorB), as.double(cache$eps.trend.alphaB), as.double(cache$eps.trend.betaB), as.double(cache$eps.beta.factorB), as.double(cache$eps.slice.factorB), as.double(cache$x_turnover0), as.double(cache$x_eps0), as.double(cache$x_turnover1), as.double(cache$x_eps1), as.double(cache$x_turnoverA), as.double(cache$x_epsA), as.double(cache$x_turnoverB), as.double(cache$x_epsB), as.double(cache$q01), as.double(cache$q10), as.double(cache$q0A), as.double(cache$qA0), as.double(cache$q1B), as.double(cache$qB1), as.double(cache$q0B), as.double(cache$qB0), as.double(cache$q1A), as.double(cache$qA1), as.double(cache$qBA), as.double(cache$qAB), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(3))
			compD[focal,] <- c(v * c(lambda0[[1]], lambda1[[1]],lambdaA[[1]], lambdaB[[1]]))
			compE[focal,] <- phi[1:4]
			if(!is.null(node)){
				fixer = c(0,0,0,0)
				fixer[state] = 1
				if(node == focal){
					compD[focal,] <- compD[focal,] * fixer
				}
				#compE[focal,] <- compE[focal,] * fixer
			}
		}else{
			lambda0 <- .C("set_birth_bisse_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$timeslice), as.double(cache$turnover.trend.alpha0), as.double(cache$turnover.trend.beta0), as.double(cache$turnover.beta.factor0), as.double(cache$turnover.slice.factor0), as.double(cache$eps.trend.alpha0), as.double(cache$eps.trend.beta0), as.double(cache$eps.beta.factor0), as.double(cache$eps.slice.factor0), as.double(cache$turnover.trend.alpha1), as.double(cache$turnover.trend.beta1), as.double(cache$turnover.beta.factor1), as.double(cache$turnover.slice.factor1), as.double(cache$eps.trend.alpha1), as.double(cache$eps.trend.beta1), as.double(cache$eps.beta.factor1), as.double(cache$eps.slice.factor1), as.double(cache$x_turnover0), as.double(cache$x_eps0), as.double(cache$x_turnover1), as.double(cache$x_eps1), as.double(cache$q01), as.double(cache$q10), as.double(cache$focal.edge.length), as.double(cache$tipward.age), as.integer(0))
			lambda1 <- .C("set_birth_bisse_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$timeslice), as.double(cache$turnover.trend.alpha0), as.double(cache$turnover.trend.beta0), as.double(cache$turnover.beta.factor0), as.double(cache$turnover.slice.factor0), as.double(cache$eps.trend.alpha0), as.double(cache$eps.trend.beta0), as.double(cache$eps.beta.factor0), as.double(cache$eps.slice.factor0), as.double(cache$turnover.trend.alpha1), as.double(cache$turnover.trend.beta1), as.double(cache$turnover.beta.factor1), as.double(cache$turnover.slice.factor1), as.double(cache$eps.trend.alpha1), as.double(cache$eps.trend.beta1), as.double(cache$eps.beta.factor1), as.double(cache$eps.slice.factor1), as.double(cache$x_turnover0), as.double(cache$x_eps0), as.double(cache$x_turnover1), as.double(cache$x_eps1), as.double(cache$q01), as.double(cache$q10), as.double(cache$focal.edge.length), as.double(cache$tipward.age), as.integer(1))
			compD[focal,] <- c(v * c(lambda0[[1]], lambda1[[1]]))
			compE[focal,] <- phi[1:2]
			if(!is.null(node)){
				if(node == focal){
					fixer = c(0,0)
					fixer[state] = 1
					compD[focal,] <- compD[focal,] * fixer
				}
			}			
		}
		###########################
		#Logcompensation bit for dealing with underflow issues. Need to give a necessary shoutout to Rich FitzJohn -- we follow his diversitree approach. VERIFIED that it works properly:
		tmp <- sum(compD[focal,])
		compD[focal,] <- compD[focal,] / tmp
		logcomp <- c(logcomp, log(tmp))
	}

	root.node <- nb.tip + 1L
	if (is.na(sum(log(compD[root.node,]))) || is.na(log(sum(1-compE[root.node,])))){
		return(bad.likelihood)
	}else{
		if(root.type == "madfitz"){
			if(hidden.states == FALSE){
				root.p = c(compD[root.node,1]/sum(compD[root.node,]), compD[root.node,2]/sum(compD[root.node,]))
			}else{
				root.p = c(compD[root.node,1]/sum(compD[root.node,]), compD[root.node,2]/sum(compD[root.node,]), compD[root.node,3]/sum(compD[root.node,]), compD[root.node,4]/sum(compD[root.node,]))
				root.p[which(is.na(root.p))] = 0 
			}
		}
		if(root.type == "equal"){
			root.p = c(rep(1/length(which(compD[root.node,] > 0)), length(compD[root.node,])))
			root.p[which(!compD[root.node,] > 0)] = 0
		}
		if(root.type == "user"){
			root.p = root.p
		}
		if(condition.on.survival == TRUE){
			if(hidden.states == FALSE){
				compD[root.node,] <- compD[root.node,] / sum(root.p * c(lambda0[[1]], lambda1[[1]]) * (1 - compE[root.node,])^2)
				#Corrects for possibility that you have 0/0:
				compD[root.node,which(is.na(compD[root.node,]))] = 0
			}else{
				compD[root.node,] <- compD[root.node,] / sum(root.p * c(lambda0[[1]], lambda1[[1]], lambdaA[[1]], lambdaB[[1]]) * (1 - compE[root.node,])^2)				
				#Corrects for possibility that you have 0/0:
				compD[root.node,which(is.na(compD[root.node,]))] = 0
			}
		}
		loglik <- log(sum(compD[root.node,] * root.p)) + sum(logcomp)
	}
	if(get.phi==TRUE){
		obj = NULL
		obj$compD.root = compD[root.node,]
		obj$compE = compE
		return(obj)
	}else{
		return(loglik)
	}
}


######################################################################################################################################
######################################################################################################################################
### Function for ascertainment bias filter -- WORK IN PROGRESS
######################################################################################################################################
######################################################################################################################################

#Rationale:
#Ascertainment bias is a common issue in biological systems. For example, for transcriptomes, shorter seqs might be easier to detect (Gao et al. 2011). Here, the clades
#selected for use in diversification analyses are biased. The bias of having to exist is already part of the model. However, we actually use more stringent filters: no one
#does a diversification analysis of a clade of Gingko: there is one species, so there is no point. The clades selected for analysis are larger than clades
#evolving for the same time with the same birth and death rates, and so rate estimates are biased towards greater net diversification. As a first attempt to deal with this
#we have allowed an ascertainment filter to be implemented. We assume that the true parameters of evolution make the examined clade exceptionally diverse in some way, and
#thus only allow parameter values that make it exceptional enough. For example, by default we assume that clades of the observed number of taxa or greater should have a 
#probability of being used of 5% or less. There is room for future improvements in dealing with this, but the effect size of this bias is large enough that it must be
#attempted to be dealt with. By changing max.probability to 1.0, the traditional approach that ignores the reality of ascertainment bias may be used.
#One possibility for this probability may be the number of taxa in the focal group divided by the number of taxa in the overall larger set of "similar" things
#For example, are whales diverse mammals? There are other mammal groups (platypus, hippos) you have chosen not to analyze, so you have to take into account your choice
#to do a moderately diverse, fairly young group. If the parameters pass this filter this, a TRUE is passed.

PassAscertainmentFilter <- function(max.probability=NULL, k, birth.rate, death.rate, time, comparison.clade.diversity=NULL, comparison.clade.age=NULL) {
	if(is.null(max.probability) && !is.null(comparison.clade.diversity)) {
		comparison.clade.rate = log(comparison.clade.diversity / 2) / comparison.clade.age
		#Magallon and Sanderson 2001 -- Eq. 4 -- we want to calculate how many lineages we expect to be around at the age of origin of the focal clade. We have only looked at our special one:
		comparison.clade.expected.lineages = 2 * exp(comparison.clade.rate * (comparison.clade.age-time)) 
		#Out of all the lineages alive at time t, we chose this one, presumably because it is diverse:
		max.probability = 1/comparison.clade.expected.lineages
	}
	net.diver.rate <- birth.rate - death.rate
	#Magallon and Sanderson 2001 -- Eq. 2a:
	exprt <- exp(net.diver.rate * time)
	beta <- (exprt - 1) / (exprt - death.rate/birth.rate)
	#Magallon and Sanderson 2001 -- Eq. 2b:
	alpha <- (death.rate/birth.rate) * beta
	#Magallon and Sanderson 2001 -- Eq. 11a:
	probNgeNtax <- ((beta^(k-2)) * ((net.diver.rate * (1 - alpha - beta + (alpha*beta)) + alpha + (2*beta) - 1))) / (1+alpha)
	return(ifelse((probNgeNtax <= max.probability) , TRUE, FALSE))
}


######################################################################################################################################
######################################################################################################################################
### Cache object for storing parameters that are used throughout hisse:
######################################################################################################################################
######################################################################################################################################

ParametersToPass <- function(phy, data, f, model.vec, timeslice, hidden.states){
	#Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPass):
	obj <- NULL
	obj$phy = phy

	if(hidden.states == FALSE){
		states = matrix(0,Ntip(phy),2)
		for(i in 1:Ntip(phy)){
			if(data[i]==0){states[i,1]=1}
			if(data[i]==1){states[i,2]=1}
			if(data[i]==2){states[i,1:2]=1}
		}
	}
	if(hidden.states == TRUE){
		states = matrix(0,Ntip(phy),4)
		for(i in 1:Ntip(phy)){
			if(data[i]==0){states[i,c(1,3)]=1}
			if(data[i]==1){states[i,c(2,4)]=1}
			if(data[i]==2){states[i,1:4]=1}
		}
	}
	if(hidden.states == "TEST"){
		states = matrix(0,Ntip(phy),4)
		for(i in 1:Ntip(phy)){
			if(data[i]==1){states[i,1]=1}
			if(data[i]==2){states[i,2]=1}
			if(data[i]==3){states[i,3]=1}
			if(data[i]==4){states[i,4]=1}
		}
	}
	obj$states = states
	obj$tot_time = max(branching.times(phy))
	obj$f = f
	if(is.null(timeslice)){
		obj$timeslice = 0
	}else{
		obj$timeslice = timeslice
	}
	obj$x_turnover0 = model.vec[1]
	obj$x_turnover1 = model.vec[2]
	obj$x_turnoverA = model.vec[3]
	obj$x_turnoverB = model.vec[4]
	
	obj$x_eps0 = model.vec[5]
	obj$x_eps1 = model.vec[6]
	obj$x_epsA = model.vec[7]
	obj$x_epsB = model.vec[8]
	
	obj$q10 = model.vec[9]
	obj$qA0 = model.vec[10]
	obj$qB0 = model.vec[11]
	obj$q01 = model.vec[12]
	obj$qA1 = model.vec[13]
	obj$qB1 = model.vec[14]
	obj$q0A = model.vec[15] 
	obj$q1A = model.vec[16]
	obj$qBA = model.vec[17] 
	obj$q0B = model.vec[18]
	obj$q1B = model.vec[19]
	obj$qAB = model.vec[20]

	obj$turnover.trend.alpha0 = model.vec[21] 
	obj$turnover.trend.alpha1 = model.vec[22] 
	obj$turnover.trend.alphaA = model.vec[23] 
	obj$turnover.trend.alphaB = model.vec[24] 
	
	obj$turnover.trend.beta0 = model.vec[25] 
	obj$turnover.trend.beta1 = model.vec[26] 
	obj$turnover.trend.betaA = model.vec[27] 
	obj$turnover.trend.betaB = model.vec[28] 
	
	obj$eps.trend.alpha0 = model.vec[29]
	obj$eps.trend.alpha1 = model.vec[30]
	obj$eps.trend.alphaA = model.vec[31]
	obj$eps.trend.alphaB = model.vec[32]

	obj$eps.trend.beta0 = model.vec[33]
	obj$eps.trend.beta1 = model.vec[34]
	obj$eps.trend.betaA = model.vec[35]
	obj$eps.trend.betaB = model.vec[36]
	
	obj$turnover.slice.factor0 = model.vec[37]
	obj$turnover.slice.factor1 = model.vec[38]
	obj$turnover.slice.factorA = model.vec[39]
	obj$turnover.slice.factorB = model.vec[40]
	
	obj$eps.slice.factor0 = model.vec[41]
	obj$eps.slice.factor1 = model.vec[42]
	obj$eps.slice.factorA = model.vec[43]
	obj$eps.slice.factorB = model.vec[44]
	
	obj$q10_slice.factor = model.vec[45]
	obj$qA0_slice.factor = model.vec[46]
	obj$qB0_slice.factor = model.vec[47]
	obj$q01_slice.factor = model.vec[48]
	obj$qA1_slice.factor = model.vec[49]
	obj$qB1_slice.factor = model.vec[50]
	obj$q0A_slice.factor = model.vec[51] 
	obj$q1A_slice.factor = model.vec[52]
	obj$qBA_slice.factor = model.vec[53] 
	obj$q0B_slice.factor = model.vec[54]
	obj$q1B_slice.factor = model.vec[55]
	obj$qAB_slice.factor = model.vec[56]
	
	obj$split.times = sort(branching.times(phy), decreasing=TRUE)
	
	return(obj)
}

######################################################################################################################################
######################################################################################################################################
### Print function for our diversity class:
######################################################################################################################################
######################################################################################################################################

##Work on this more
print.hisse.fit <- function(x,...){
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips,row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
	cat("Probability of extinction:", x$phi, "\n")
	cat("\n")
	param.est0 <- matrix(0,4,2)
	param.est0[,1] <- x$solution[c(1,21,25,37)]
	param.est0[,2] <- x$solution[c(5,29,33,41)]
	if(x$hidden.states == TRUE){
		param.est0 <- data.frame(param.est0, row.names=c("rate0A", "alpha0A", "beta0A", "timeslice.factor0A"))
	}else{
		param.est0 <- data.frame(param.est0, row.names=c("rate0", "alpha0", "beta0", "timeslice.factor0"))
		
	}
	param.est1 <- matrix(0,4,2)
	param.est1[,1] <- x$solution[c(2,22,26,38)]
	param.est1[,2] <- x$solution[c(6,30,34,42)]
	if(x$hidden.states == TRUE){
		param.est1 <- data.frame(param.est1, row.names=c("rate1A", "alpha1A", "beta1A", "timeslice.factor1A"))
	}else{
		param.est1 <- data.frame(param.est1, row.names=c("rate1", "alpha1", "beta1", "timeslice.factor1"))
	}
	
	names(param.est0) <- names(param.est1) <- c("turnover", "extinction")	
	
	param.est.sp.0 <- param.est0[1,1] / (1 + param.est0[1,2])
	param.est.mu.0 <- (param.est0[1,1] * param.est0[1,2]) / (1 + param.est0[1,2])
	param.est.sp.1 <- param.est1[1,1] / (1 + param.est1[1,2])
	param.est.mu.1 <- (param.est1[1,1] * param.est1[1,2]) / (1 + param.est1[1,2])

	if(x$output.type == "net.div"){
		param.est0[1,1] <- param.est.sp.0 - param.est.mu.0
		param.est1[1,1] <- param.est.sp.1 - param.est.mu.1
		names(param.est0) <- names(param.est1) <- c("net.div", "extinction")	
	}
	
	if(x$output.type == "raw"){
		param.est0[1,1:2] <- c(param.est.sp.0, param.est.mu.0)
		param.est1[1,1:2] <- c(param.est.sp.1, param.est.mu.1)
		names(param.est0) <- names(param.est1) <- c("lambda", "mu")	
	}
	
	cat("Diversification Rates\n")
	print(param.est0)
	cat("\n")
	print(param.est1)
	cat("\n")
	if(x$hidden.states == TRUE){
		param.estA <- matrix(0,4,2)
		param.estA[,1] <- x$solution[c(3,23,27,39)]
		param.estA[,2] <- x$solution[c(7,31,35,43)]
		param.estA <- data.frame(param.estA, row.names=c("rate0B", "alpha0B", "beta0B", "timeslice.factor0B"))
		param.estB <- matrix(0,4,2)
		param.estB[,1] <- x$solution[c(4,24,28,49)]
		param.estB[,2] <- x$solution[c(8,32,36,44)]
		param.estB <- data.frame(param.estB, row.names=c("rate1B", "alpha1B", "beta1B", "timeslice.factor1B"))
		names(param.estA) <- names(param.estB) <- c("turnover", "extinction")	
		
		param.est.sp.A <- param.estA[1,1] / (1 + param.estA[1,2])
		param.est.mu.A <- (param.estA[1,1] * param.estA[1,2]) / (1 + param.estA[1,2])
		param.est.sp.B <- param.estB[1,1] / (1 + param.estB[1,2])
		param.est.mu.B <- (param.estB[1,1] * param.estB[1,2]) / (1 + param.estB[1,2])

		if(x$output.type == "net.div"){
			param.estA[1,1] <- param.est.sp.A - param.est.mu.A
			param.estB[1,1] <- param.est.sp.B - param.est.mu.B
			names(param.estA) <- names(param.estB) <- c("net.div", "extinction")	
		}

		if(x$output.type == "raw"){
			param.estA[1,1:2] <- c(param.est.sp.A, param.est.mu.A)
			param.estB[1,1:2] <- c(param.est.sp.B, param.est.mu.B)
			names(param.estA) <- names(param.estB) <- c("lambda", "mu")	
		}
		
		print(param.estA)
		cat("\n")
		print(param.estB)
		cat("\n")
		index.mat<-matrix(NA,4,4)
		diag(index.mat)<-13
		index.mat[is.na(index.mat)]<-1:12
		t.rates <- x$solution[c(9:20)]
		q.mat <- matrix(t.rates[index.mat],dim(index.mat))
		rownames(q.mat) <- c("(0A)","(1A)","(0B)","(1B)")
		colnames(q.mat) <- c("(0A)","(1A)","(0B)","(1B)")			
		cat("Transition Rates\n")
		print(q.mat)
		cat("\n")
	}else{
		q.mat <- matrix(0,2,2)
		q.mat[2,1] <- x$solution[9]
		q.mat[1,2] <- x$solution[12]
		rownames(q.mat) <- c("(0)","(1)")
		colnames(q.mat) <- c("(0)","(1)")					
		cat("Transition Rates\n")
		print(q.mat)
		cat("\n")
	}
}

