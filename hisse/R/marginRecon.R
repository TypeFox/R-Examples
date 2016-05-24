######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION -- calculates marginal probabilities
######################################################################################################################################
######################################################################################################################################

MarginRecon <- function(phy, data, f, pars, hidden.states=TRUE, four.state.null=FALSE, timeslice=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, verbose=TRUE, n.cores=NULL){
	if(!is.null(root.p)) {
		root.type="user"
		root.p <- root.p / sum(root.p)	
		if(hidden.states ==TRUE & length(root.p)==2){
			root.p <- rep(root.p, 2)
			root.p <- root.p / sum(root.p)	
			warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
		}
	}
	
	if(four.state.null == FALSE){
		phy$node.label = NULL
		data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
		data.new <- data.new[phy$tip.label,]
		model.vec = pars
		#Prerequisites for running the downpass algorithm:
		cache = ParametersToPass(phy, data.new[,1], model.vec, f=f, timeslice=NULL, hidden.states=hidden.states)
		cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
		cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
		cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
		cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])
		cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
		cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
		cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
		cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])	
		
		nb.tip <- length(phy$tip.label)
		nb.node <- phy$Nnode
		if(hidden.states == FALSE){
			nstates=2
		}else{
			nstates=4
		}
		nodes <- unique(phy$edge[,1])	
		
		if(is.null(n.cores)){
			marginal.probs <- matrix(0, nb.node+nb.tip, nstates)
			for (i in seq(from = 1, length.out = nb.node)) {
				focal <- nodes[i]
				marginal.probs.tmp <- c()
				for (j in 1:nstates){
					marginal.probs.tmp <- c(marginal.probs.tmp, DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
				}
				best.probs = max(marginal.probs.tmp)
				marginal.probs.rescaled = marginal.probs.tmp - best.probs
				marginal.probs[focal,] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
				if (verbose && i%%100==0) {
					print(paste(i, "of", nb.node, "nodes done"))	
				}
			}
			if(hidden.states==TRUE){
				for (i in seq(from = 1, length.out = nb.tip)) {
					marginal.probs.tmp <- numeric(4)
					nstates = which(!cache$states[i,] == 0)
					for (j in nstates){
						cache$states.keep = cache$states[i,]
						cache$states[i,] = 0
						cache$states[i,j] = 1
						marginal.probs.tmp[j] <- DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
					}
					cache$states[i,] = cache$states.keep
					best.probs = max(marginal.probs.tmp[nstates])
					marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
					marginal.probs[i,nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
					if (verbose && i%%100==0) {
						print(paste(i, "of", nb.tip, "tips done"))	
					}			
				}
			}
			obj <- NULL
			if(!hidden.states == TRUE){
				marginal.probs[1:nb.tip,] = cache$states
			}
			marginal.probs <- cbind(1:(nb.node+nb.tip), marginal.probs)
			obj$node.mat = marginal.probs[-(1:nb.tip),] 
			obj$tip.mat = marginal.probs[1:nb.tip,]
			if(hidden.states == TRUE){
				raw.rates <- ParameterTransform(x=model.vec[1:4], y=model.vec[5:8])
				rates.mat <- matrix(0, 5, 4)
				rates.mat[1,] <- model.vec[1:4]
				rates.mat[2,] <- raw.rates[1:4] - raw.rates[5:8]
				rates.mat[3,] <- raw.rates[1:4]
				rates.mat[4,] <- model.vec[5:8]
				rates.mat[5,] <- raw.rates[5:8]
				rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
				colnames(rates.mat) <- c("0A", "1A", "0B", "1B")
				colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "1A", "0B", "1B")
			}else{
				raw.rates <- ParameterTransform(x=model.vec[1:4], y=model.vec[5:8])
				rates.mat <- matrix(0, 5, 2)
				rates.mat[1,] <- model.vec[1:2]
				rates.mat[2,] <- raw.rates[1:2] - raw.rates[5:6]
				rates.mat[3,] <- raw.rates[1:2]
				rates.mat[4,] <- model.vec[5:6]
				rates.mat[5,] <- raw.rates[5:6]
				rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
				colnames(rates.mat) <- c("0", "1")
				colnames(obj$node.mat) <- colnames(obj$tip.mat) <- c("id", "0", "1")
			}
			obj$rates.mat = rates.mat			
			phy$node.label = apply(marginal.probs[,-1], 1, which.max)[-(1:nb.tip)]
			obj$phy = phy
		}else{
			NodeEval <- function(node){
				focal <- node
				marginal.probs.tmp <- c()
				for (j in 1:nstates){
					marginal.probs.tmp <- c(marginal.probs.tmp, DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
				}
				best.probs = max(marginal.probs.tmp)
				marginal.probs.rescaled = marginal.probs.tmp - best.probs
				marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
				return(c(node, marginal.probs))
			}
			node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
			
			if(hidden.states==TRUE){
				TipEval <- function(tip){
					marginal.probs.tmp <- numeric(4)
					nstates = which(!cache$states[tip,] == 0)
					for (j in nstates){
						cache$states.keep = cache$states[tip,]
						cache$states[tip,] = 0
						cache$states[tip,j] = 1
						marginal.probs.tmp[j] <- DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
					}
					cache$states[tip,] = cache$states.keep
					best.probs = max(marginal.probs.tmp[nstates])
					marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
					marginal.probs <- numeric(4)
					marginal.probs[nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
					return(c(tip, marginal.probs))
				}
				tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
			}
			obj <- NULL
			if(hidden.states == TRUE){
				obj$node.mat <- matrix(unlist(node.marginals), ncol = 5, byrow = TRUE)
				obj$tip.mat = matrix(unlist(tip.marginals), ncol = 5, byrow = TRUE)
				raw.rates <- ParameterTransform(x=model.vec[1:4], y=model.vec[5:8])
				rates.mat <- matrix(0, 5, 4)
				rates.mat[1,] <- model.vec[1:4]
				rates.mat[2,] <- raw.rates[1:4] - raw.rates[5:8]
				rates.mat[3,] <- raw.rates[1:4]
				rates.mat[4,] <- model.vec[5:8]
				rates.mat[5,] <- raw.rates[5:8]
				rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
				colnames(rates.mat) <- c("0A", "1A", "0B", "1B")
				colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "1A", "0B", "1B")
			}else{
				obj$node.mat <- matrix(unlist(node.marginals), ncol = 3, byrow = TRUE)
				obj$tip.mat = cbind(1:Ntip(phy), cache$states)
				raw.rates <- ParameterTransform(x=model.vec[1:4], y=model.vec[5:8])
				rates.mat <- matrix(0, 5, 2)
				rates.mat[1,] <- model.vec[1:2]
				rates.mat[2,] <- raw.rates[1:2] - raw.rates[5:6]
				rates.mat[3,] <- raw.rates[1:2]
				rates.mat[4,] <- model.vec[5:6]
				rates.mat[5,] <- raw.rates[5:6]
				rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
				colnames(rates.mat) <- c("0", "1")
				colnames(obj$node.mat) <- colnames(obj$tip.mat) <- c("id", "0", "1")
			}
			obj$rates.mat = rates.mat
			phy$node.label = apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
			obj$phy = phy
		}
	}else{
		phy$node.label = NULL
		data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
		data.new <- data.new[phy$tip.label,]
		model.vec = pars
		#Prerequisites for running the downpass algorithm:
		cache = ParametersToPassNull(phy, data.new[,1], model.vec, f=f)
		
		nb.tip <- length(phy$tip.label)
		nb.node <- phy$Nnode
		nstates=8
		nodes <- unique(phy$edge[,1])	
		
		if(is.null(n.cores)){
			marginal.probs <- matrix(0, nb.node+nb.tip, nstates)
			for (i in seq(from = 1, length.out = nb.node)) {
				focal <- nodes[i]
				marginal.probs.tmp <- c()
				for (j in 1:nstates){
					marginal.probs.tmp <- c(marginal.probs.tmp, DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
				}
				best.probs = max(marginal.probs.tmp)
				marginal.probs.rescaled = marginal.probs.tmp - best.probs
				marginal.probs[focal,] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
				if (verbose && i%%100==0) {
					print(paste(i, "of", nb.node, "nodes done"))	
				}
			}
			#Now for the tips...
			for (i in seq(from = 1, length.out = nb.tip)) {
				marginal.probs.tmp <- numeric(4)
				nstates = which(!cache$states[i,] == 0)
				for (j in nstates){
					cache$states.keep = cache$states[i,]
					cache$states[i,] = 0
					cache$states[i,j] = 1
					marginal.probs.tmp[j] <- DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
				}
				cache$states[i,] = cache$states.keep
				best.probs = max(marginal.probs.tmp[nstates])
				marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
				marginal.probs[i,nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
				if (verbose && i%%100==0) {
					print(paste(i, "of", nb.tip, "tips done"))	
				}			
			}
			obj <- NULL
			obj$node.mat = marginal.probs[-(1:nb.tip),] 
			obj$tip.mat = marginal.probs[1:nb.tip,]
			raw.rates <- ParameterTransform(x=model.vec[1:8], y=model.vec[9:16])
			rates.mat <- matrix(0, 5, 8)
			rates.mat[1,] <- model.vec[1:8]
			rates.mat[2,] <- raw.rates[1:8] - raw.rates[9:16]
			rates.mat[3,] <- raw.rates[1:8]
			rates.mat[4,] <- model.vec[9:16]
			rates.mat[5,] <- raw.rates[9:16]
			rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
			colnames(rates.mat) <- c("0A", "0B", "0C", "0D","1A", "1B", "1C", "1D")
			colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "0B", "0C", "0D","1A", "1B", "1C", "1D")
			obj$rates.mat = rates.mat			
			phy$node.label = apply(marginal.probs, 1, which.max)[-(1:nb.tip)]
			obj$phy = phy
		}else{
			NodeEval <- function(node){
				focal <- node
				marginal.probs.tmp <- c()
				for (j in 1:nstates){
					marginal.probs.tmp <- c(marginal.probs.tmp, DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
				}
				best.probs = max(marginal.probs.tmp)
				marginal.probs.rescaled = marginal.probs.tmp - best.probs
				marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
				return(c(node, marginal.probs))
			}
			node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
			
			TipEval <- function(tip){
				marginal.probs.tmp <- numeric(8)
				nstates = which(!cache$states[tip,] == 0)
				for (j in nstates){
					cache$states.keep = cache$states[tip,]
					cache$states[tip,] = 0
					cache$states[tip,j] = 1
					marginal.probs.tmp[j] <- DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
				}
				cache$states[tip,] = cache$states.keep
				best.probs = max(marginal.probs.tmp[nstates])
				marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
				marginal.probs <- numeric(8)
				marginal.probs[nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
				return(c(tip, marginal.probs))
			}
			tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
			obj <- NULL
			obj$node.mat <- matrix(unlist(node.marginals), ncol = 9, byrow = TRUE)
			obj$tip.mat = matrix(unlist(tip.marginals), ncol = 9, byrow = TRUE)
			raw.rates <- ParameterTransform(x=model.vec[1:8], y=model.vec[9:16])
			rates.mat <- matrix(0, 5, 8)
			rates.mat[1,] <- model.vec[1:8]
			rates.mat[2,] <- raw.rates[1:8] - raw.rates[9:16]
			rates.mat[3,] <- raw.rates[1:8]
			rates.mat[4,] <- model.vec[9:16]
			rates.mat[5,] <- raw.rates[9:16]
			rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
			colnames(rates.mat) <- c("0A", "0B", "0C", "0D","1A", "1B", "1C", "1D")
			colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "0B", "0C", "0D","1A", "1B", "1C", "1D")
			obj$rates.mat = rates.mat
			phy$node.label = apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
			obj$phy = phy
		}
	}
	if(!is.null(aic)){
		obj$aic = aic
	}
	class(obj) = "hisse.states"
	return(obj)
}


ParameterTransform <- function(x, y){
	speciation <- x / (1+y)
	extinction <- (x*y) / (1+y)
	return(c(speciation, extinction))
}


print.hisse.states <- function(x,...){
	print(x$phy)
}



