#HIDDEN RATES MODEL OF BINARY TRAIT EVOLUTION

#written by Jeremy M. Beaulieu

corHMM<-function(phy, data, rate.cat, rate.mat=NULL, node.states=c("joint", "marginal", "scaled"), optim.method=c("subplex"), p=NULL, root.p=NULL, ip=NULL, nstarts=10, n.cores=NULL, sann.its=5000, diagn=FALSE){
	
	# Checks to make sure node.states is not NULL.  If it is, just returns a diagnostic message asking for value.
	if(is.null(node.states)){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to corHMM command for parameter \'node.states\': joint, marginal, or scaled.")
		return(obj)
	}
	else { # even if node.states is not NULL, need to make sure its one of the three valid options
		valid.models <- c("joint", "marginal", "scaled")
		if(!any(valid.models == node.states)){
			obj <- NULL
			obj$loglik <- NULL
			obj$diagnostic <- paste("\'",node.states, "\' is not valid for ancestral state reconstruction method.  Please pass one of the following to corHMM command for parameter \'node.states\': joint, marginal, or scaled.",sep="")
			return(obj)
		}
		if(length(node.states) > 1){ # User did not enter a value, so just pick marginal.
			node.states <- "marginal"
			cat("No model selected for \'node.states\'. Will perform marginal ancestral state estimation.\n")
		}
	}
	
	# Checks to make sure phy & data have same taxa. Fixes conflicts (see match.tree.data function).
	matching <- match.tree.data(phy,data) 
	data <- matching$data
	phy <- matching$phy
	
	# Will not perform reconstructions on invariant characters
	if(nlevels(as.factor(data[,1])) <= 1){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("Character is invariant. Analysis stopped.",sep="")
		return(obj)
	} else {
		# Still need to make sure second level isnt just an ambiguity
		lvls <- as.factor(data[,1])
		if(nlevels(as.factor(data[,1])) == 2 && length(which(lvls == "?"))){
			obj <- NULL
			obj$loglik <- NULL
			obj$diagnostic <- paste("Character is invariant. Analysis stopped.",sep="")
			return(obj)
		}
	}
	
	#Creates the data structure and orders the rows to match the tree. 
	phy$edge.length[phy$edge.length<=1e-5]=1e-5
	data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
	data.sort <- data.sort[phy$tip.label,]
	
	counts <- table(data.sort[,1])
	levels <- levels(as.factor(data.sort[,1]))
	cols <- as.factor(data.sort[,1])
	cat("State distribution in data:\n")
	cat("States:",levels,"\n",sep="\t")
	cat("Counts:",counts,"\n",sep="\t")
	
	#Some initial values for use later
	k=2
	
	ub = log(100)
	lb = -20
	
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	rate.cat=rate.cat
	root.p=root.p
	nstarts=nstarts
	ip=ip
	
	model.set.final<-rate.cat.set.corHMM(phy=phy,data.sort=data.sort,rate.cat=rate.cat)
	if(!is.null(rate.mat)){
		rate <- rate.mat
		model.set.final$np <- max(rate, na.rm=TRUE)
		rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
		model.set.final$rate <- rate
		model.set.final$index.matrix <- rate.mat
	}

	lower = rep(lb, model.set.final$np)
	upper = rep(ub, model.set.final$np)
	
	if(optim.method=="twoStep"){
		if(!is.null(p)){
			cat("Calculating likelihood from a set of fixed parameters", "\n")
			out<-NULL
			est.pars<-log(p)
			out$objective <- dev.corhmm(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			est.pars <- exp(est.pars)
		}
		else{
			cat("Initializing...", "\n")
			model.set.init<-rate.cat.set.corHMM(phy=phy,data.sort=data.sort,rate.cat=1)
			rate<-rate.mat.maker(hrm=TRUE,rate.cat=1)
			rate<-rate.par.eq(rate,eq.par=c(1,2))
			model.set.init$index.matrix<-rate
			rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			model.set.init$rate<-rate
			dat<-as.matrix(data.sort)
			dat<-phyDat(dat,type="USER", levels=c("0","1"))
			par.score<-parsimony(phy, dat, method="fitch")
			tl <- sum(phy$edge.length)
			mean.change = par.score/tl
			if(mean.change==0){
				ip=exp(lb)+0.01
			}else{
				ip<-rexp(1, 1/mean.change)
			}
			if(ip < exp(lb) || ip > exp(ub)){ # initial parameter value is outside bounds
				ip <- exp(lb)
			}			
			lower = rep(lb, model.set.final$np)
			upper = rep(ub, model.set.final$np)
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
			cat("Finished. Beginning simulated annealing Round 1...", "\n")
			out.sann <- GenSA(rep(log(ip), model.set.final$np), fn=dev.corhmm, lower=lower, upper=upper, control=list(max.call=sann.its), phy=phy,liks=model.set.final$liks, Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			cat("Finished. Refining using subplex routine...", "\n")
            #out = subplex(out.sann$par, fn=dev.corhmm, control=list(.Machine$double.eps^0.25, parscale=rep(0.1, length(out.sann$par))), phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p)
			out = nloptr(x0=out.sann$par, eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			cat("Finished. Beginning simulated annealing Round 2...", "\n")
			out.sann <- GenSA(out$solution, fn=dev.corhmm, lower=lower, upper=upper, control=list(max.call=sann.its), phy=phy,liks=model.set.final$liks, Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			cat("Finished. Refining using subplex routine...", "\n")
            #out = subplex(out.sann$par, fn=dev.corhmm, control=list(.Machine$double.eps^0.25, parscale=rep(0.1, length(out.sann$par))), phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p)
			out = nloptr(x0=out.sann$par, eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			cat("Finished. Beginning simulated annealing Round 3...", "\n")
			out.sann <- GenSA(out$solution, fn=dev.corhmm, lower=lower, upper=upper, control=list(max.call=sann.its), phy=phy,liks=model.set.final$liks, Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			cat("Finished. Refining using subplex routine...", "\n")
            #out = subplex(out.sann$par, fn=dev.corhmm, control=list(.Machine$double.eps^0.25, parscale=rep(0.1, length(out.sann$par))), phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p)
			out = nloptr(x0=out.sann$par, eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)

			loglik <- -out$objective
			est.pars <- exp(out$solution)
		}
	}
	if(optim.method=="subplex"){
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
		if(!is.null(p)){
			cat("Calculating likelihood from a set of fixed parameters", "\n")
			out<-NULL
			est.pars<-log(p)
			out$objective<-dev.corhmm(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			loglik <- -out$objective
			est.pars <- exp(est.pars)
		}
		#If a user-specified starting value(s) is not supplied this begins loop through a set of randomly chosen starting values:
		else{
			#If a user-specified starting value(s) is supplied:
			if(is.null(ip)){
				if(is.null(n.cores)){
					cat("Beginning thorough optimization search -- performing", nstarts, "random restarts", "\n")
					#If the analysis is to be run a single processor:
					if(is.null(n.cores)){
						#Sets parameter settings for random restarts by taking the parsimony score and dividing
						#by the total length of the tree
						dat<-as.matrix(data.sort)
						dat<-phyDat(dat,type="USER", levels=c("0","1"))
						par.score<-parsimony(phy, dat, method="fitch")/2
						tl <- sum(phy$edge.length)
						mean.change = par.score/tl
						if(mean.change==0){
							ip=0.01+exp(lb)
						}else{
							ip <- rexp(model.set.final$np, 1/mean.change)
						}
						ip[ip < exp(lb)] = exp(lb)
						ip[ip > exp(ub)] = exp(lb)
						out = nloptr(x0=rep(log(ip), length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
						tmp = matrix(,1,ncol=(1+model.set.final$np))
						tmp[,1] = out$objective
						tmp[,2:(model.set.final$np+1)] = out$solution
						for(i in 2:nstarts){
						#Temporary solution for ensuring an ordered Q with respect to the rate classes. If a simpler model is called this feature is automatically turned off:
							if(mean.change==0){
								starts=runif(0.01+exp(lb), 1, model.set.final$np)
							}else{
								starts<-rexp(model.set.final$np, 1/mean.change)
							}
							starts[starts < exp(lb)] = exp(lb)
							starts[starts > exp(ub)] = exp(lb)							
							par.order<-NA
							if(rate.cat == 2){
								try(par.order<-starts[3] > starts[8])
								if(!is.na(par.order)){
									pp.tmp <- c(starts[3],starts[8])
									starts[3] <- min(pp.tmp)
									starts[8] <- max(pp.tmp)
								}
							}
							if(rate.cat == 3){
								try(par.order <- starts[3] > starts[9] | starts[9] > starts[14])
								if(!is.na(par.order)){
									pp.tmp <- c(starts[3],starts[9],starts[14])
									starts[3] <- min(pp.tmp)
									starts[9] <- median(pp.tmp)
									starts[14] <- max(pp.tmp)
								}
							}
							if(rate.cat == 4){
								try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[20])
								if(!is.na(par.order)){
									pp.tmp <- c(starts[3],starts[9],starts[15],starts[20])
									starts[3] <- pp.tmp[order(pp.tmp)][1]
									starts[9] <- pp.tmp[order(pp.tmp)][2]
									starts[15] <- pp.tmp[order(pp.tmp)][3]
									starts[20] <- pp.tmp[order(pp.tmp)][4]
								}
							}
							if(rate.cat == 5){
								try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[21] | starts[21] > starts[26])
								if(!is.na(par.order)){
									pp.tmp <- c(starts[3],starts[9],starts[15],starts[21],starts[26])
									starts[3] <- pp.tmp[order(pp.tmp)][1]
									starts[9] <- pp.tmp[order(pp.tmp)][2]
									starts[15] <- pp.tmp[order(pp.tmp)][3]
									starts[21] <- pp.tmp[order(pp.tmp)][4]
									starts[26] <- pp.tmp[order(pp.tmp)][5]	
								}
							}
							if(rate.cat == 6){
								try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[21] | starts[21] > starts[27] | starts[27] > starts[32])
								if(!is.na(par.order)){
									pp.tmp <- c(starts[3],starts[9],starts[15],starts[21],starts[27],starts[32])
									starts[3] <- pp.tmp[order(pp.tmp)][1]
									starts[9] <- pp.tmp[order(pp.tmp)][2]
									starts[15] <- pp.tmp[order(pp.tmp)][3]
									starts[21] <- pp.tmp[order(pp.tmp)][4]
									starts[27] <- pp.tmp[order(pp.tmp)][5]
									starts[32] <- pp.tmp[order(pp.tmp)][6]
								}
							}
							if(rate.cat == 7){
								try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[21] | starts[21] > starts[27] | starts[27] > starts[33] | starts[33] > starts[38])
								if(!is.na(par.order)){
									pp.tmp <- c(starts[3],starts[9],starts[15],starts[21],starts[27],starts[33],starts[38])
									starts[3] <- pp.tmp[order(pp.tmp)][1]
									starts[9] <- pp.tmp[order(pp.tmp)][2]
									starts[15] <- pp.tmp[order(pp.tmp)][3]
									starts[21] <- pp.tmp[order(pp.tmp)][4]
									starts[27] <- pp.tmp[order(pp.tmp)][5]
									starts[33] <- pp.tmp[order(pp.tmp)][6]
									starts[38] <- pp.tmp[order(pp.tmp)][7]
								}
							}						
							out.alt = nloptr(x0=rep(log(starts), length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
							tmp[,1] = out.alt$objective
							tmp[,2:(model.set.final$np+1)] = starts
							if(out.alt$objective < out$objective){
								out = out.alt
								ip = starts
							}
							else{
								out = out
								ip = ip
							}
						}
						loglik <- -out$objective
						est.pars <- exp(out$solution)
					}
				}
				#If the analysis is to be run on multiple processors:
				else{
					#Sets parameter settings for random restarts by taking the parsimony score and dividing
					#by the total length of the tree
					cat("Beginning thorough optimization search -- performing", nstarts, "random restarts", "\n")
					dat<-as.matrix(data.sort)
					dat<-phyDat(dat,type="USER", levels=c("0","1"))
					par.score<-parsimony(phy, dat, method="fitch")/2
					tl <- sum(phy$edge.length)
					mean.change = par.score/tl
					random.restart<-function(nstarts){
						tmp = matrix(,1,ncol=(1+model.set.final$np))
						#Temporary solution for ensuring an ordered Q with respect to the rate classes. If a simpler model is called this feature is automatically turned off:
						if(mean.change==0){
							starts=rep(0.01+exp(lb), model.set.final$np)
						}else{
							starts<-rexp(model.set.final$np, 1/mean.change)
						}
						starts[starts < exp(lb)] = exp(lb)
						starts[starts > exp(ub)] = exp(lb)
						par.order<-NA
						if(rate.cat == 2){
							try(par.order<-starts[3] > starts[8])
							if(!is.na(par.order)){
								pp.tmp <- c(starts[3],starts[8])
								starts[3] <- min(pp.tmp)
								starts[8] <- max(pp.tmp)
							}
						}
						if(rate.cat == 3){
							try(par.order <- starts[3] > starts[9] | starts[9] > starts[14])
							if(!is.na(par.order)){
								pp.tmp <- c(starts[3],starts[9],starts[14])
								starts[3] <- min(pp.tmp)
								starts[9] <- median(pp.tmp)
								starts[14] <- max(pp.tmp)
							}
						}
						if(rate.cat == 4){
							try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[20])
							if(!is.na(par.order)){
								pp.tmp <- c(starts[3],starts[9],starts[15],starts[20])
								starts[3] <- pp.tmp[order(pp.tmp)][1]
								starts[9] <- pp.tmp[order(pp.tmp)][2]
								starts[15] <- pp.tmp[order(pp.tmp)][3]
								starts[20] <- pp.tmp[order(pp.tmp)][4]
							}
						}
						if(rate.cat == 5){
							try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[21] | starts[21] > starts[26])
							if(!is.na(par.order)){
								pp.tmp <- c(starts[3],starts[9],starts[15],starts[21],starts[26])
								starts[3] <- pp.tmp[order(pp.tmp)][1]
								starts[9] <- pp.tmp[order(pp.tmp)][2]
								starts[15] <- pp.tmp[order(pp.tmp)][3]
								starts[21] <- pp.tmp[order(pp.tmp)][4]
								starts[26] <- pp.tmp[order(pp.tmp)][5]	
							}
						}
						if(rate.cat == 6){
							try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[21] | starts[21] > starts[27] | starts[27] > starts[32])
							if(!is.na(par.order)){
								pp.tmp <- c(starts[3],starts[9],starts[15],starts[21],starts[27],starts[32])
								starts[3] <- pp.tmp[order(pp.tmp)][1]
								starts[9] <- pp.tmp[order(pp.tmp)][2]
								starts[15] <- pp.tmp[order(pp.tmp)][3]
								starts[21] <- pp.tmp[order(pp.tmp)][4]
								starts[27] <- pp.tmp[order(pp.tmp)][5]
								starts[32] <- pp.tmp[order(pp.tmp)][6]
							}
						}
						if(rate.cat == 7){
							try(par.order <- starts[3] > starts[9] | starts[9] > starts[15] | starts[15] > starts[21] | starts[21] > starts[27] | starts[27] > starts[33] | starts[33] > starts[38])
							if(!is.na(par.order)){
								pp.tmp <- c(starts[3],starts[9],starts[15],starts[21],starts[27],starts[33],starts[38])
								starts[3] <- pp.tmp[order(pp.tmp)][1]
								starts[9] <- pp.tmp[order(pp.tmp)][2]
								starts[15] <- pp.tmp[order(pp.tmp)][3]
								starts[21] <- pp.tmp[order(pp.tmp)][4]
								starts[27] <- pp.tmp[order(pp.tmp)][5]
								starts[33] <- pp.tmp[order(pp.tmp)][6]
								starts[38] <- pp.tmp[order(pp.tmp)][7]
							}
						}
						out = nloptr(x0=log(starts), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy, liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
						tmp[,1] = out$objective
						tmp[,2:(model.set.final$np+1)] = out$solution
						tmp
					}
					restart.set<-mclapply(1:nstarts, random.restart, mc.cores=n.cores)
					#Finds the best fit within the restart.set list
					best.fit<-which.min(unlist(lapply(1:nstarts,function(i) lapply(restart.set[[i]][,1],min))))
					#Generates an object to store results from restart algorithm:
					out<-NULL
					out$objective=unlist(restart.set[[best.fit]][,1])
					out$solution=unlist(restart.set[[best.fit]][,2:(model.set.final$np+1)])
					loglik <- -out$objective
					est.pars <- exp(out$solution)
				}
			}
			else{
				cat("Beginning subplex optimization routine -- Starting value(s):", ip, "\n")
				ip=ip
				out = nloptr(x0=rep(log(ip), length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
				loglik <- -out$objective
				est.pars <- exp(out$solution)
			}			
		}
	}
	
	#Starts the summarization process:
	cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
	
	TIPS <- 1:nb.tip
	if (node.states == "marginal" || node.states == "scaled"){
		lik.anc <- ancRECON(phy, data, est.pars, hrm=TRUE, rate.cat, rate.mat=rate.mat, method=node.states, ntraits=NULL, root.p=root.p)
		pr<-apply(lik.anc$lik.anc.states,1,which.max)
		phy$node.label <- pr
		tip.states <- lik.anc$lik.tip.states
		row.names(tip.states) <- phy$tip.label
	}
	if (node.states == "joint"){
		lik.anc <- ancRECON(phy, data, est.pars, hrm=TRUE, rate.cat, rate.mat=rate.mat, method=node.states, ntraits=NULL,root.p=root.p)
		phy$node.label <- lik.anc$lik.anc.states
		tip.states <- lik.anc$lik.tip.states
	}
	
	cat("Finished. Performing diagnostic tests.", "\n")
	
	#Approximates the Hessian using the numDeriv function
	if(diagn==TRUE){
		h <- hessian(func=dev.corhmm, x=est.pars, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
		solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
		solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix], dim(model.set.final$index.matrix))
		hess.eig <- eigen(h,symmetric=TRUE)
		eigval<-signif(hess.eig$values, 2)
		eigvect<-round(hess.eig$vectors, 2)		
	}
	else{
		solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
		solution.se <- matrix(0,dim(solution)[1],dim(solution)[1])
		eigval<-NULL
		eigvect<-NULL
	}
	if (rate.cat == 1){
		rownames(solution) <- rownames(solution.se) <- c("(0)","(1)")
		colnames(solution) <- colnames(solution.se) <- c("(0)","(1)")			
		#Initiates user-specified reconstruction method:
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){
				colnames(lik.anc$lik.anc.states) <- c("P(0)","P(1)")
			}
		}
	}
	if (rate.cat == 2){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){		
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
			}
		}
	}
	if (rate.cat == 3){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){		
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
			}
		}
	}
	if (rate.cat == 4){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){	
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
			}
		}
	}
	if (rate.cat == 5){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){	
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
			}
		}
	}
	if (rate.cat == 6){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)")
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){	
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)")
			}
		}
	}
	if (rate.cat == 7){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)","(0,R7)","(1,R7)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)","(0,R7)","(1,R7)")
		if (is.character(node.states)) {
			if (node.states == "marginal" || node.states == "scaled"){	
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)","(0,R7)","(1,R7)")
			}
		}
	}
	obj = list(loglik = loglik, AIC = -2*loglik+2*model.set.final$np,AICc = -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))),rate.cat=rate.cat,solution=solution, solution.se=solution.se, index.mat=model.set.final$index.matrix, data=data.sort, phy=phy, states=lik.anc$lik.anc.states, tip.states=tip.states, iterations=out$iterations, eigval=eigval, eigvect=eigvect) 
	class(obj)<-"corhmm"
	return(obj)
}

#Print function
print.corhmm<-function(x,...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","Rate.cat","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	param.est<- x$solution
	cat("Rates\n")
	print(param.est)
	cat("\n")
	
	if(any(x$eigval<0)){
		index.matrix <- x$index.mat
		#If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
		if (any(x$eigval<0)) {
			cat("The objective function may be at a saddle point", "\n")
		}
	}
	else{
		cat("Arrived at a reliable solution","\n")
	}
}

#Generalized ace() function that allows analysis to be carried out when there are polytomies:
dev.corhmm <- function(p,phy,liks,Q,rate,root.p) {
	p = exp(p)
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	TIPS <- 1:nb.tip
	comp <- numeric(nb.tip + nb.node)
	phy <- reorder(phy, "pruningwise")
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])
	k.rates <- dim(Q)[2] / 2
	
	if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
	
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	for (i  in seq(from = 1, length.out = nb.node)) {
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		v <- 1
		#Loops through all descendants of focal (how we deal with polytomies):
		for (desIndex in sequence(length(desRows))){
			v<-v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
		}
		#Sum the likelihoods:
		comp[focal] <- sum(v)
		#Divide each likelihood by the sum to obtain probabilities:
		liks[focal, ] <- v/comp[focal]
	}
	#Temporary solution for ensuring an ordered Q with respect to the rate classes. If a simpler model is called this feature is automatically turned off:
	par.order<-NA
	if(k.rates == 2){
		try(par.order <- p[3] > p[8])
		if(!is.na(par.order)){
			if(par.order == TRUE){
				return(1000000)
			}
		}
	}
	if(k.rates == 3){
		try(par.order <- p[3] > p[9] | p[9] > p[14])
		if(!is.na(par.order)){
			if(par.order == TRUE){
				return(1000000)
			}
		}		
	}
	if(k.rates == 4){
		try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[20])
		if(!is.na(par.order)){
			if(par.order == TRUE){
				return(1000000)
			}
		}		
	}
	if(k.rates == 5){
		try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[26])
		if(!is.na(par.order)){
			if(par.order == TRUE){
				return(1000000)
			}
		}		
	}
	if(k.rates == 6){
		try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[27] | p[27] > p[32])
		if(!is.na(par.order)){
			if(par.order == TRUE){
				return(1000000)
			}
		}		
	}
	if(k.rates == 7){
		try(par.order <- p[3] > p[9] | p[9] > p[15] | p[15] > p[21] | p[21] > p[27] | p[27] > p[33] | p[33] > p[38])
		if(!is.na(par.order)){
			if(par.order == TRUE){
				return(1000000)
			}
		}		
	}	
	#Specifies the root:
	root <- nb.tip + 1L
	#If any of the logs have NAs restart search:
	if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
	else{
		equil.root <- NULL
		for(i in 1:ncol(Q)){
			posrows <- which(Q[,i] >= 0)
			rowsum <- sum(Q[posrows,i])
			poscols <- which(Q[i,] >= 0)
			colsum <- sum(Q[i,poscols])
			equil.root <- c(equil.root,rowsum/(rowsum+colsum))
		}		
		if (is.null(root.p)){
            flat.root = equil.root
			k.rates <- 1/length(which(!is.na(equil.root)))
			flat.root[!is.na(flat.root)] = k.rates
			flat.root[is.na(flat.root)] = 0
			loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))
		}
		else{
			if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)
                if(root.p == "yang"){
                    diag(Q) = 0
                    equil.root = colSums(Q) / sum(Q)
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(equil.root)+log(liks[root,])))))
                    if(is.infinite(loglik)){
                        return(1000000)
                    }
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks[root,] / sum(liks[root,])
                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
                }
			}
			# root.p!==NULL will fix root probabilities based on user supplied vector:
			else{
                loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))
				if(is.infinite(loglik)){
                    return(1000000)
                }
			}
		}
	}
	loglik
}

rate.cat.set.corHMM<-function(phy,data.sort,rate.cat){
	
	k=2
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	obj$rate.cat<-rate.cat
	
	rate<-rate.mat.maker(hrm=TRUE,rate.cat=rate.cat)
	index.matrix<-rate
	rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
	
	#Makes a matrix of tip states and empty cells corresponding 
	#to ancestral nodes during the optimization process.	
	x <- data.sort[,1]
	TIPS <- 1:nb.tip
	
	for(i in 1:nb.tip){
		if(is.na(x[i])){x[i]=2}
	}
	if (rate.cat == 1){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,1]=1}
			if(x[i]==1){liks[i,2]=1}
			if(x[i]==2){liks[i,1:2]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 2){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3)]=1}
			if(x[i]==1){liks[i,c(2,4)]=1}
			if(x[i]==2){liks[i,1:4]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 3){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5)]=1}
			if(x[i]==1){liks[i,c(2,4,6)]=1}
			if(x[i]==2){liks[i,1:6]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 4){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8)]=1}
			if(x[i]==2){liks[i,1:8]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 5){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7,9)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8,10)]=1}
			if(x[i]==2){liks[i,1:10]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 6){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7,9,11)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8,10,12)]=1}
			if(x[i]==2){liks[i,1:12]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 7){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7,9,11,13)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8,10,12,14)]=1}
			if(x[i]==2){liks[i,1:14]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	obj$np<-max(rate)-1
	obj$rate<-rate
	obj$index.matrix<-index.matrix
	obj$liks<-liks
	obj$Q<-Q
	
	obj
}

