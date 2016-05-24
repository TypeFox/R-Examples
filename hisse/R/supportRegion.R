
######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap -- Simulating confidence intervals for parameters estimated in HiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegion <- function(hisse.obj, n.points=1000, scale.int=0.1, desired.delta=2, output.type="turnover", hidden.states=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, verbose=TRUE){
	phy <- hisse.obj$phy
	data <- hisse.obj$data
	data.new<-data.frame(data[,2], data[,2], row.names=data[,1])
	data.new<-data.new[phy$tip.label,]
	f <- hisse.obj$f
	np <- max(hisse.obj$index.par)-1
	par <- numeric(np)
	free.parameters <- which(hisse.obj$index.par < max(hisse.obj$index.par))
	np.sequence <- 1:np
	for(i in np.sequence){
		par[i] <- hisse.obj$solution[which(hisse.obj$index.par == np.sequence[i])[1]]
	}
	lower <- rep(exp(-20), np)
	upper <- rep(exp(20), np)

	#Bad Jeremy! Hard-coded column headers...
    if(output.type == "turnover"){
        interval.names <- c("lnLik", "turn.0A", "turn.1A", "turn.0B", "turn.1B", "eps.0A", "eps.1A", "eps.0B", "eps.1B","q1A0A","q0B0A","q1B0A","q0A1A","q0B1A","q1B1A","q0A0B","q1A0B","q1B0B","q0A1B","q1A1B","q0A1B","turn.alpha.0A","turn.alpha.1A", "turn.alpha.0B", "turn.alpha.1B", "turn.beta.0A","turn.beta.1A", "turn.beta.0B", "turn.beta.1B", "eps.alpha.0A","eps.alpha.1A", "eps.alpha.0B", "eps.alpha.1B", "eps.beta.0A","eps.beta.1A", "eps.beta.0B", "eps.beta.1B", "turn.slice.0A","turn.slice.1A", "turn.slice.0B", "turn.slice.1B", "eps.slice.0A","eps.slice.1A", "eps.slice.0B", "eps.slice.1B", "q0A1A.slice","q1A0A.slice","q0A0B.slice","q0B0A.slice","q1A1B.slice","q1B1A.slice","q0A1B.slice","q1B0A.slice","q1A0B.slice","q0B1A.slice","q1B0B.slice","q0B1B.slice")
    }
    if(output.type == "net.div"){
        interval.names <- c("lnLik", "netdiv.0A", "netdiv.1A", "netdiv.0B", "netdiv.1B", "eps.0A", "eps.1A", "eps.0B", "eps.1B","q1A0A","q0B0A","q1B0A","q0A1A","q0B1A","q1B1A","q0A0B","q1A0B","q1B0B","q0A1B","q1A1B","q0A1B","turn.alpha.0A","turn.alpha.1A", "turn.alpha.0B", "turn.alpha.1B", "turn.beta.0A","turn.beta.1A", "turn.beta.0B", "turn.beta.1B", "eps.alpha.0A","eps.alpha.1A", "eps.alpha.0B", "eps.alpha.1B", "eps.beta.0A","eps.beta.1A", "eps.beta.0B", "eps.beta.1B", "turn.slice.0A","turn.slice.1A", "turn.slice.0B", "turn.slice.1B", "eps.slice.0A","eps.slice.1A", "eps.slice.0B", "eps.slice.1B", "q0A1A.slice","q1A0A.slice","q0A0B.slice","q0B0A.slice","q1A1B.slice","q1B1A.slice","q0A1B.slice","q1B0A.slice","q1A0B.slice","q0B1A.slice","q1B0B.slice","q0B1B.slice")
    }
    if(output.type == "raw"){
        interval.names <- c("lnLik", "lambda.0A", "lambda.1A", "lambda.0B", "lambda.1B", "mu.0A", "mu.1A", "mu.0B", "mu.1B","q1A0A","q0B0A","q1B0A","q0A1A","q0B1A","q1B1A","q0A0B","q1A0B","q1B0B","q0A1B","q1A1B","q0A1B","turn.alpha.0A","turn.alpha.1A", "turn.alpha.0B", "turn.alpha.1B", "turn.beta.0A","turn.beta.1A", "turn.beta.0B", "turn.beta.1B", "eps.alpha.0A","eps.alpha.1A", "eps.alpha.0B", "eps.alpha.1B", "eps.beta.0A","eps.beta.1A", "eps.beta.0B", "eps.beta.1B", "turn.slice.0A","turn.slice.1A", "turn.slice.0B", "turn.slice.1B", "eps.slice.0A","eps.slice.1A", "eps.slice.0B", "eps.slice.1B", "q0A1A.slice","q1A0A.slice","q0A0B.slice","q0B0A.slice","q1A1B.slice","q1B1A.slice","q0A1B.slice","q1B0A.slice","q1A0B.slice","q0B1A.slice","q1B0B.slice","q0B1B.slice")
    }

	interval.results <- AdaptiveConfidenceIntervalSampling(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=hisse.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int)
    interval.results.final <- matrix(0, n.points+1, length(hisse.obj$index.par))
	for(i in 1:(n.points+1)){
		par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
		interval.results.final[i,] <- c(par.rep,0)[hisse.obj$index.par]
	}
	interval.results.final[,21:56] = 1
	interval.results.final <- cbind(interval.results[,1], interval.results.final)
    if(output.type == "net.div"){
        lambda.0A <- interval.results.final[,2] / (1 + interval.results.final[,6])
        lambda.1A <- interval.results.final[,3] / (1 + interval.results.final[,7])
        lambda.0B <- interval.results.final[,4] / (1 + interval.results.final[,8])
        lambda.1B <- interval.results.final[,5] / (1 + interval.results.final[,9])
        mu.0A <- (interval.results.final[,2] * interval.results.final[,6]) / (1 + interval.results.final[,6])
        mu.1A <- (interval.results.final[,3] * interval.results.final[,7]) / (1 + interval.results.final[,7])
        mu.0B <- (interval.results.final[,4] * interval.results.final[,8]) / (1 + interval.results.final[,8])
        mu.1B <- (interval.results.final[,5] * interval.results.final[,9]) / (1 + interval.results.final[,9])
        interval.results.final[,2] <- lambda.0A - mu.0A
        interval.results.final[,3] <- lambda.1A - mu.1A
        interval.results.final[,4] <- lambda.0B - mu.0B
        interval.results.final[,5] <- lambda.1B - mu.1B
    }
    if(output.type == "raw"){
        lambda.0A <- interval.results.final[,2] / (1 + interval.results.final[,6])
        lambda.1A <- interval.results.final[,3] / (1 + interval.results.final[,7])
        lambda.0B <- interval.results.final[,4] / (1 + interval.results.final[,8])
        lambda.1B <- interval.results.final[,5] / (1 + interval.results.final[,9])
        mu.0A <- (interval.results.final[,2] * interval.results.final[,6]) / (1 + interval.results.final[,6])
        mu.1A <- (interval.results.final[,3] * interval.results.final[,7]) / (1 + interval.results.final[,7])
        mu.0B <- (interval.results.final[,4] * interval.results.final[,8]) / (1 + interval.results.final[,8])
        mu.1B <- (interval.results.final[,5] * interval.results.final[,9]) / (1 + interval.results.final[,9])
        interval.results.final[,2] <- lambda.0A
        interval.results.final[,3] <- lambda.1A
        interval.results.final[,4] <- lambda.0B
        interval.results.final[,5] <- lambda.1B
        interval.results.final[,6] <- mu.0A
        interval.results.final[,7] <- mu.1A
        interval.results.final[,8] <- mu.0B
        interval.results.final[,9] <- mu.1B
    }
	interval.results.in <- interval.results.final[which(interval.results.final[,1] - min(interval.results.final[,1])<=desired.delta),]
    ci.interval = apply(interval.results.in, 2, quantile)
	colnames(interval.results.final) <- colnames(interval.results.in) <- colnames(ci.interval) <- interval.names

	obj = NULL
	obj$ci <- ci.interval[,1:21]
	obj$points.within.region = interval.results.in[,1:21]
	obj$all.points = interval.results.final[,1:21]
	class(obj) = "hisse.support"	
	return(obj)
}


AdaptiveConfidenceIntervalSampling <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, index.par, f, hidden.states, condition.on.survival, root.type, root.p, scale.int) {
	
	#Wrangle the data so that we can make use of DownPass easily:
	actual.params = which(index.par < max(index.par))
	model.vec <- numeric(length(index.par))
	model.vec[] <- c(par,0)[index.par]
	model.vec.tmp = model.vec[21:56]
	model.vec.tmp[model.vec.tmp==0] = 1
	model.vec[21:56] = model.vec.tmp		
	cache = ParametersToPass(phy, data[,1], model.vec, f=f, timeslice=NULL, hidden.states=hidden.states) 
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
	cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
	cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])	
	cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
	cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
	cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
	cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])
	phy$node.label <- NULL
	#############################################################
	#Now assess the likelihood at the MLE:
	starting <- -DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
    #Generate the multipliers for feeling the boundaries:
	min.multipliers <- rep(1, length(par))
	max.multipliers <- rep(1, length(par))
	results <- data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
	results[1,] <- unname(c(starting, par))
	for (i in sequence(n.points)) {
		sim.points <- NA
		while(is.na(sim.points[1])) {
			sim.points <- GenerateValues(par, lower=lower, upper=upper, scale.int=scale.int, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
		}
		par <- sim.points
		model.vec <- numeric(length(index.par))
		model.vec[] <- c(sim.points,0)[index.par]
		model.vec.tmp = model.vec[21:56]
		model.vec.tmp[model.vec.tmp==0] = 1
		model.vec[21:56] = model.vec.tmp		
		cache = ParametersToPass(phy, data[,1], model.vec, f=f, timeslice=NULL, hidden.states=hidden.states) 
		cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
		cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
		cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
		cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])	
		cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
		cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
		cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
		cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])
		second <- -DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
		results[i+1,] <- c(second, sim.points)
		if(i%%20==0) {
			for (j in sequence(length(par))) {
				returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
				total.range <- range(results[,j+1], na.rm=TRUE)
				width.ratio <- diff(returned.range)/diff(total.range)
				if(is.na(width.ratio)) {
					width.ratio=1	
				}
				if(width.ratio > 0.5) { #we are not sampling widely enough
					min.multipliers[j] <- min.multipliers[j] * (1-scale.int)
					max.multipliers[j] <- max.multipliers[j] * (1+scale.int) #expand the range
				} else {
					min.multipliers[j] <- 1
					max.multipliers[j] <- 1
				}
			}
		}
		if (verbose && i%%100==0) {
			print(paste(i, "of", n.points, "points done"))	
		}
	}
	return(results)
}


GenerateValues <- function(par, lower, upper, scale.int, max.tries=100, expand.prob=0, examined.max, examined.min) {
	pass=FALSE
	tries=0
	while(!pass && tries<=max.tries) {
		tries <- tries+1
		pass=TRUE
		new.vals <- rep(NA, length(par))
		for(i in sequence(length(par))) {
			examined.max[i] <- max(0.001, examined.max[i])
			new.vals[i] <- runif(1, max(lower[i], (1-scale.int)*examined.min[i]), min(upper[i], (1+scale.int)*examined.max[i]))
			if(new.vals[i]<lower[i]) {
				pass=FALSE
			}
			if(new.vals[i]>upper[i]) {
				pass=FALSE
			}
		}
	}
	if(tries>max.tries) {
		return(NA)
	}
	return(new.vals)
}


print.hisse.support <- function(x,...){

	cat("\nSupport Region\n")
	print(x$ci[,-1])
	cat("\n")
}
	


