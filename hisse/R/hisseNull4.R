
######################################################################################################################################
######################################################################################################################################
### HiSSE four state null -- A null HiSSE model that assumes four possible hidden states, none associated with a character state
######################################################################################################################################
######################################################################################################################################

hisse.null4 <- function(phy, data, f=c(1,1), turnover.anc=rep(c(1,2,3,4),2), eps.anc=rep(c(1,2,3,4),2), trans.type = "equal", condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, output.type="turnover", sann=FALSE, sann.its=10000, max.tol=.Machine$double.eps^.25){
	
	#Some basic formatting of parameters:
	phy$node.label <- NULL	
	sub.mat1 <- sub.mat2 <- TransMatMaker(hidden.states=TRUE)
	sub.mat3 <- sub.mat4 <- matrix(NA, 4,4)
	if(trans.type == "equal"){
		diag(sub.mat3) <- diag(sub.mat4) <- 1
		trans.mat <-rbind(cbind(sub.mat1, sub.mat3),cbind(sub.mat4,sub.mat2))
		trans.mat[!is.na(trans.mat)] = 1
		rownames(trans.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")
		colnames(trans.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")			
	}else{
		sub.mat1 <- TransMatMaker(hidden.states=TRUE)
		sub.mat2 <- TransMatMaker(hidden.states=TRUE)
		sub.mat1[!is.na(sub.mat1)] <- sub.mat2[!is.na(sub.mat2)] <- 1
		diag(sub.mat4) <- 2
		diag(sub.mat3) <- 3
		trans.mat <-rbind(cbind(sub.mat1,sub.mat3),cbind(sub.mat4,sub.mat2))
		rownames(trans.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")
		colnames(trans.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")			
	}
	trans.rate <- trans.mat

	#Initialize our pars vector and build from there:
	pars = c(turnover.anc)
	#Add in extinction fraction:
	eps.anc.tmp = eps.anc
	eps.anc.tmp[which(eps.anc.tmp>0)] = (eps.anc.tmp[which(eps.anc.tmp>0)] + max(pars))
	pars = c(pars, eps.anc.tmp)
	#Add in transition rates:
	trans.rate.tmp = trans.rate[!is.na(trans.rate)]
	trans.rate.tmp[which(trans.rate.tmp > 0)] = (trans.rate.tmp[which(trans.rate.tmp > 0)] + max(pars))
	pars = c(pars, trans.rate.tmp)
	
	np <- max(pars)
	pars[pars==0] <- np+1
	
	cat("Initializing...", "\n")
	
	data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
	data.new <- data.new[phy$tip.label,]

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
		def.set.pars <- c(rep(log(init.pars[1]+init.pars[3]), 8), rep(log(init.pars[3]/init.pars[1]),8), rep(log(init.pars[5]), 32))
	}else{
		init.pars <- starting.point.generator(phy, 2, samp.freq.tree, yule=FALSE)
		names(init.pars) <- NULL
		init.eps = init.pars[3]/init.pars[1]
		if(init.eps == 0){
			init.eps = 1e-6
		}
		def.set.pars <- c(rep(log(init.pars[1]+init.pars[3]), 8), rep(log(init.eps),8), rep(log(init.pars[5]), 32))
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
		out = subplex(ip, fn=DevOptimizeNull, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, phy=phy, data=data.new[,1], f=f, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np)
	}else{
		cat("Finished. Beginning simulated annealing...", "\n")
		out.sann = GenSA(ip, fn=DevOptimizeNull, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, phy=phy, data=data.new[,1], f=f, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np)
		cat("Finished. Refining using subplex routine...", "\n")
		out = subplex(out.sann$par, fn=DevOptimizeNull, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, phy=phy, data=data.new[,1], f=f, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np)
	}
	solution <- numeric(length(pars))
	solution[] <- c(exp(out$par), 0)[pars]
	loglik = -out$value
	
	cat("Finished. Summarizing results...", "\n")
	
	obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy, data=data, output.type=output.type, trans.type=trans.type, trans.mat=trans.mat, max.tol=max.tol) 
	class(obj) = "hisse.null4.fit"		
	
	return(obj)		
}


######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################

#Function used for optimizing parameters:
DevOptimizeNull <- function(p, pars, phy, data, f, condition.on.survival, root.type, root.p, timeslice, np) {
	#Generates the final vector with the appropriate parameter estimates in the right place:
	p.new <- exp(p)
	model.vec <- numeric(length(pars))
	model.vec[] <- c(p.new, 0)[pars]
	cache = ParametersToPassNull(phy, data, model.vec, f=f)
	logl <- DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
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

DownPassNull <- function(phy, cache, bad.likelihood=-10000000000, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL) {
	#Some preliminaries:
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	phy <- reorder(phy, "pruningwise")
	anc <- unique(phy$edge[,1])		
	TIPS <- 1:nb.tip

	compD <- matrix(0, nrow=nb.tip + nb.node, ncol=8)
	compE <- matrix(0, nrow=nb.tip + nb.node, ncol=8)		
		
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
		v = rep(1, 8)
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
			pars <- list(cache$tot_time, cache$timeslice, cache$x_turnover0A, cache$x_turnover0B, cache$x_turnover0C, cache$x_turnover0D, cache$x_turnover1A, cache$x_turnover1B, cache$x_turnover1C, cache$x_turnover1D, cache$x_eps0A, cache$x_eps0B, cache$x_eps0C, cache$x_eps0D, cache$x_eps1A, cache$x_eps1B, cache$x_eps1C, cache$x_eps1D, cache$q0B0A, cache$q0C0A, cache$q0D0A, cache$q1A0A, cache$q0A0B, cache$q0C0B, cache$q0D0B, cache$q1B0B, cache$q0A0C, cache$q0B0C, cache$q0D0C, cache$q1C0C, cache$q0A0D, cache$q0B0D, cache$q0C0D, cache$q1D0D, cache$q0A1A, cache$q1B1A, cache$q1C1A, cache$q1D1A, cache$q0B1B, cache$q1A1B, cache$q1C1B, cache$q1D1B, cache$q0C1C, cache$q1A1C, cache$q1B1C, cache$q1D1C, cache$q0D1D, cache$q1A1D, cache$q1B1D, cache$q1C1D, cache$focal_edge_length, cache$tipward_age)
			NUMELEMENTS <- 51 #needed for passing in vector to C
			padded.pars <- rep(0, NUMELEMENTS)
			pars <- c(unlist(pars))
			stopifnot(length(padded.pars)<=NUMELEMENTS)
			padded.pars[sequence(length(pars))]<-pars				
			yini <-c(E0A=cache$node.E[1], E0B=cache$node.E[2], E0C=cache$node.E[3], E0D=cache$node.E[4], E1A=cache$node.E[5], E1B=cache$node.E[6], E1C=cache$node.E[7], E1D=cache$node.E[8], D0A=cache$node.D[1], D0B=cache$node.D[2], D0C=cache$node.D[3], D0D=cache$node.D[4], D1A=cache$node.D[5], D1B=cache$node.D[6], D1C=cache$node.D[7], D1D=cache$node.D[8])
			times=c(cache$tipward.age, cache$rootward.age)
			prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_hisse_null", padded.pars, initfunc="initmod_hisse_null", dllname = "hisse", rtol=1e-8, atol=1e-8, hmax=10)

			######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
			if(attributes(prob.subtree.cal.full)$istate[1] < 0){
				return(bad.likelihood)
			}else{
				prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
			}
			##############################################################################
			if(!is.finite(prob.subtree.cal[9]) | !is.finite(prob.subtree.cal[10]) | !is.finite(prob.subtree.cal[11]) | !is.finite(prob.subtree.cal[12]) | !is.finite(prob.subtree.cal[13]) | !is.finite(prob.subtree.cal[14]) | !is.finite(prob.subtree.cal[15]) | !is.finite(prob.subtree.cal[16])){
				return(bad.likelihood)
			}				
			if(prob.subtree.cal[9]<0 | prob.subtree.cal[10]<0 | prob.subtree.cal[11]<0 | prob.subtree.cal[12]<0 | prob.subtree.cal[13]<0 | prob.subtree.cal[14]<0 | prob.subtree.cal[15]<0 | prob.subtree.cal[16]<0){
					return(bad.likelihood)
			}
			#Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
			phi <- c(phi, prob.subtree.cal[1:8])
			v <- v * prob.subtree.cal[9:16]				
		}
		#C call to set_birth_void -- NOTE: The first input is zero as we need to declare the birth_rate. It gets written over and is now the first element in the list that is returned. Everything else should be self explanatory.
		lambda0A <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(0))
		lambda0B <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(1))
		lambda0C <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(2))
		lambda0D <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(3))
		lambda1A <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(4))
		lambda1B <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(5))
		lambda1C <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(6))
		lambda1D <- .C("set_birth_hisse_null_void", as.double(0.0), as.double(cache$rootward.age), as.double(cache$tot_time), as.double(cache$x_turnover0A), as.double(cache$x_eps0A), as.double(cache$x_turnover0B), as.double(cache$x_eps0B),as.double(cache$x_turnover0C), as.double(cache$x_eps0C), as.double(cache$x_turnover0D), as.double(cache$x_eps0D), as.double(cache$x_turnover1A), as.double(cache$x_eps1A), as.double(cache$x_turnover1B), as.double(cache$x_eps1B),as.double(cache$x_turnover1C), as.double(cache$x_eps1C), as.double(cache$x_turnover1D), as.double(cache$x_eps1D), as.double(cache$q0B0A), as.double(cache$q0C0A), as.double(cache$q0D0A), as.double(cache$q1A0A), as.double(cache$q0A0B), as.double(cache$q0C0B), as.double(cache$q0D0B), as.double(cache$q1B0B), as.double(cache$q0A0C), as.double(cache$q0B0C), as.double(cache$q0D0C), as.double(cache$q1C0C), as.double(cache$q0A0D), as.double(cache$q0B0D), as.double(cache$q0C0D), as.double(cache$q1D0D), as.double(cache$q0A1A), as.double(cache$q1B1A), as.double(cache$q1C1A), as.double(cache$q1D1A), as.double(cache$q0B1B), as.double(cache$q1A1B), as.double(cache$q1C1B), as.double(cache$q1D1B), as.double(cache$q0C1C), as.double(cache$q1A1C), as.double(cache$q1B1C), as.double(cache$q1D1C), as.double(cache$q0D1D), as.double(cache$q1A1D), as.double(cache$q1B1D), as.double(cache$q1C1D), as.double(cache$focal.edge.ength), as.double(cache$tipward.age), as.integer(7))	
		compD[focal,] <- c(v * c(lambda0A[[1]], lambda0B[[1]],lambda0C[[1]], lambda0D[[1]], lambda1A[[1]], lambda1B[[1]], lambda1C[[1]], lambda1D[[1]]))
		compE[focal,] <- phi[1:8]
		if(!is.null(node)){
			fixer = rep(0, 8)
			fixer[state] = 1
			if(node == focal){
				compD[focal,] <- compD[focal,] * fixer
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
			root.p = c(compD[root.node,1]/sum(compD[root.node,]), compD[root.node,2]/sum(compD[root.node,]), compD[root.node,3]/sum(compD[root.node,]), compD[root.node,4]/sum(compD[root.node,]), compD[root.node,5]/sum(compD[root.node,]), compD[root.node,6]/sum(compD[root.node,]), compD[root.node,7]/sum(compD[root.node,]), compD[root.node,8]/sum(compD[root.node,]))
			root.p[which(is.na(root.p))] = 0 
		}
		if(root.type == "equal"){
			root.p = c(rep(1/length(which(compD[root.node,] > 0)), length(compD[root.node,])))
			root.p[which(!compD[root.node,] > 0)] = 0
		}
		if(root.type == "user"){
			root.p = root.p
		}
		if(condition.on.survival == TRUE){
			compD[root.node,] <- compD[root.node,] / sum(root.p * c(lambda0A[[1]], lambda0B[[1]], lambda0C[[1]], lambda0D[[1]], lambda1A[[1]], lambda1B[[1]], lambda1C[[1]], lambda1D[[1]]) * (1 - compE[root.node,])^2)		
			#Corrects for possibility that you have 0/0:
			compD[root.node,which(is.na(compD[root.node,]))] = 0
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

#PassAscertainmentFilter <- function(max.probability=NULL, k, birth.rate, death.rate, time, comparison.clade.diversity=NULL, comparison.clade.age=NULL) {
#	if(is.null(max.probability) && !is.null(comparison.clade.diversity)) {
#		comparison.clade.rate = log(comparison.clade.diversity / 2) / comparison.clade.age
#		#Magallon and Sanderson 2001 -- Eq. 4 -- we want to calculate how many lineages we expect to be around at the age of origin of the focal clade. We have only looked at our special one:
#		comparison.clade.expected.lineages = 2 * exp(comparison.clade.rate * (comparison.clade.age-time)) 
#		#Out of all the lineages alive at time t, we chose this one, presumably because it is diverse:
#		max.probability = 1/comparison.clade.expected.lineages
#	}
#	net.diver.rate <- birth.rate - death.rate
#	#Magallon and Sanderson 2001 -- Eq. 2a:
#	exprt <- exp(net.diver.rate * time)
#	beta <- (exprt - 1) / (exprt - death.rate/birth.rate)
#	#Magallon and Sanderson 2001 -- Eq. 2b:
#	alpha <- (death.rate/birth.rate) * beta
#	#Magallon and Sanderson 2001 -- Eq. 11a:
#	probNgeNtax <- ((beta^(k-2)) * ((net.diver.rate * (1 - alpha - beta + (alpha*beta)) + alpha + (2*beta) - 1))) / (1+alpha)
#	return(ifelse((probNgeNtax <= max.probability) , TRUE, FALSE))
#}


######################################################################################################################################
######################################################################################################################################
### Cache object for storing parameters that are used throughout hisse.null4:
######################################################################################################################################
######################################################################################################################################

ParametersToPassNull <- function(phy, data, f, model.vec){
	#Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassNull):
	obj <- NULL
	obj$phy = phy

	states = matrix(0,Ntip(phy),8)
	for(i in 1:Ntip(phy)){
		if(data[i]==0){states[i,c(1,2,3,4)]=1}
		if(data[i]==1){states[i,c(5,6,7,8)]=1}
		if(data[i]==2){states[i,1:8]=1}
	}
	
	obj$states = states
	obj$tot_time = max(branching.times(phy))
	obj$f = f

	obj$x_turnover0A = model.vec[1]
	obj$x_turnover0B = model.vec[2]
	obj$x_turnover0C = model.vec[3]
	obj$x_turnover0D = model.vec[4]
	obj$x_turnover1A = model.vec[5]
	obj$x_turnover1B = model.vec[6]
	obj$x_turnover1C = model.vec[7]
	obj$x_turnover1D = model.vec[8]
	
	obj$x_eps0A = model.vec[9]
	obj$x_eps0B = model.vec[10]
	obj$x_eps0C = model.vec[11]
	obj$x_eps0D = model.vec[12]
	obj$x_eps1A = model.vec[13]
	obj$x_eps1B = model.vec[14]
	obj$x_eps1C = model.vec[15]
	obj$x_eps1D = model.vec[16]

	#So many transitions...
	obj$q0B0A = model.vec[17]
	obj$q0C0A = model.vec[18]
	obj$q0D0A = model.vec[19]
	obj$q1A0A = model.vec[20]
	obj$q0A0B = model.vec[21]
	obj$q0C0B = model.vec[22]
	obj$q0D0B = model.vec[23]
	obj$q1B0B = model.vec[24]
	obj$q0A0C = model.vec[25]
	obj$q0B0C = model.vec[26] 
	obj$q0D0C = model.vec[27]
	obj$q1C0C = model.vec[28]
	obj$q0A0D = model.vec[29]
	obj$q0B0D = model.vec[30] 
	obj$q0C0D = model.vec[31]
	obj$q1D0D = model.vec[32]
	obj$q0A1A = model.vec[33]
	obj$q1B1A = model.vec[34] 
	obj$q1C1A = model.vec[35]
	obj$q1D1A = model.vec[36]
	obj$q0B1B = model.vec[37]
	obj$q1A1B = model.vec[38] 
	obj$q1C1B = model.vec[39]
	obj$q1D1B = model.vec[40]
	obj$q0C1C = model.vec[41]
	obj$q1A1C = model.vec[42] 
	obj$q1B1C = model.vec[43]
	obj$q1D1C = model.vec[44]
	obj$q0D1D = model.vec[45]
	obj$q1A1D = model.vec[46] 
	obj$q1B1D = model.vec[47]
	obj$q1C1D = model.vec[48]
	
	obj$split.times = sort(branching.times(phy), decreasing=TRUE)
	
	return(obj)
}

######################################################################################################################################
######################################################################################################################################
### Print function for our hisse.null4 class:
######################################################################################################################################
######################################################################################################################################

print.hisse.null4.fit <- function(x,...){
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips,row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
	param.est <- matrix(0,8,2)
	param.est[,1] <- x$solution[c(1:8)]
	param.est[,2] <- x$solution[c(9:16)]
	param.est <- data.frame(param.est, row.names=c("rate0A", "rate0B", "rate0C", "rate0D", "rate1A", "rate1B", "rate1C", "rate1D"))
	
	names(param.est) <- c("turnover", "extinction")	
	
	param.est.sp.0A <- param.est[1,1] / (1 + param.est[1,2])
	param.est.mu.0A <- (param.est[1,1] * param.est[1,2]) / (1 + param.est[1,2])
	param.est.sp.0B <- param.est[2,1] / (1 + param.est[2,2])
	param.est.mu.0B <- (param.est[2,1] * param.est[2,2]) / (1 + param.est[2,2])
	param.est.sp.0C <- param.est[3,1] / (1 + param.est[3,2])
	param.est.mu.0C <- (param.est[3,1] * param.est[3,2]) / (1 + param.est[3,2])
	param.est.sp.0D <- param.est[4,1] / (1 + param.est[4,2])
	param.est.mu.0D <- (param.est[4,1] * param.est[4,2]) / (1 + param.est[4,2])

	param.est.sp.1A <- param.est[5,1] / (1 + param.est[5,2])
	param.est.mu.1A <- (param.est[5,1] * param.est[5,2]) / (1 + param.est[5,2])
	param.est.sp.1B <- param.est[6,1] / (1 + param.est[6,2])
	param.est.mu.1B <- (param.est[6,1] * param.est[6,2]) / (1 + param.est[6,2])
	param.est.sp.1C <- param.est[7,1] / (1 + param.est[7,2])
	param.est.mu.1C <- (param.est[7,1] * param.est[7,2]) / (1 + param.est[7,2])
	param.est.sp.1D <- param.est[8,1] / (1 + param.est[8,2])
	param.est.mu.1D <- (param.est[8,1] * param.est[8,2]) / (1 + param.est[8,2])
	
	if(x$output.type == "net.div"){
		param.est[1,1] <- param.est.sp.0A - param.est.mu.0A
		param.est[2,1] <- param.est.sp.0B - param.est.mu.0B
		param.est[3,1] <- param.est.sp.0C - param.est.mu.0C
		param.est[4,1] <- param.est.sp.0D - param.est.mu.0D
		param.est[5,1] <- param.est.sp.1A - param.est.mu.1A
		param.est[6,1] <- param.est.sp.1B - param.est.mu.1B
		param.est[7,1] <- param.est.sp.1C - param.est.mu.1C
		param.est[8,1] <- param.est.sp.1D - param.est.mu.1D
		names(param.est) <- c("net.div", "extinction")	
	}
	
	if(x$output.type == "raw"){
		param.est[1,1:2] <- c(param.est.sp.0A, param.est.mu.0A)
		param.est[2,1:2] <- c(param.est.sp.0B, param.est.mu.0B)
		param.est[3,1:2] <- c(param.est.sp.0C, param.est.mu.0C)
		param.est[4,1:2] <- c(param.est.sp.0D, param.est.mu.0D)
		param.est[5,1:2] <- c(param.est.sp.1A, param.est.mu.1A)
		param.est[6,1:2] <- c(param.est.sp.1B, param.est.mu.1B)
		param.est[7,1:2] <- c(param.est.sp.1C, param.est.mu.1C)
		param.est[8,1:2] <- c(param.est.sp.1D, param.est.mu.1D)		
		names(param.est0) <- names(param.est1) <- c("lambda", "mu")	
	}
	
	cat("Diversification Rates\n")
	print(param.est)
	cat("\n")
	x$trans.mat[!is.na(x$trans.mat)] = 1:32
	x$trans.mat[is.na(x$trans.mat)] = 33
	t.rates <- x$solution[c(17:48)]
	q.mat <- matrix(t.rates[x$trans.mat],dim(x$trans.mat))
	rownames(q.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")
	colnames(q.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")			
	cat("Transition Rates\n")
	print(q.mat)
}

