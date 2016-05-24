#DISCRETE TRAIT TREE PAINTING METHOD##

#written by Jeremy M. Beaulieu

corPAINT<-function(phy, data, ntraits=2, rate.mat=NULL, model=c("ER","SYM","ARD"), node.states=c("joint", "marginal", "scaled"), p=NULL, root.p=NULL, ip=NULL, lb=0, ub=100, diagn=FALSE){
	
	if(is.null(phy$node.label)==TRUE){
		stop("There are no node label designations")
	}
	
	#Creates the data structure and orders the rows to match the tree
	phy$edge.length[phy$edge.length==0]=1e-5

	if(ntraits==1){
		data.sort<-data.frame(data[,2], data[,3], row.names=data[,1])
	}	
	if(ntraits==2){
		data.sort<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
	}
	if(ntraits==3){
		data.sort<-data.frame(data[,2], data[,3], data[,4], data[,5], row.names=data[,1])
	}
	data.sort<-data.sort[phy$tip.label,]
	#Some initial values for use later - will clean up
	k=ntraits
	nl=2
	
	# Check to make sure values are reasonable (i.e. non-negative)
	if(ub < 0){
		ub <- 100
	}
	if(lb < 0){
		lb <- 0
	}
	if(ub < lb){ # This user really needs help
		ub <- 100
		lb <- 0
	}

	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	
	ntraits=ntraits
	model=model
	root.p=root.p
	ip=ip

	if(ntraits==1){
		tot.states <- factor(c(phy$node.label,as.character(data[,3])))
		nregimes <- length(levels(tot.states))
		regimes <- numeric(nb.tip+nb.node)
		TIPS <- 1:nb.tip		
		regimes[TIPS]<-data[,3]
		regimes[-TIPS]<-as.numeric(phy$node.label)
	}
	if(ntraits==2){
		tot.states <- factor(c(phy$node.label,as.character(data[,4])))
		nregimes <- length(levels(tot.states))
		regimes <- numeric(nb.tip+nb.node)
		TIPS <- 1:nb.tip		
		regimes[TIPS]<-data[,4]
		regimes[-TIPS]<-as.numeric(phy$node.label)
	}
	if(ntraits==3){
		tot.states <- factor(c(phy$node.label,as.character(data[,5])))
		nregimes <- length(levels(tot.states))
		regimes <- numeric(nb.tip+nb.node)
		TIPS <- 1:nb.tip				
		regimes[TIPS] <- data[,5]
		regimes[-TIPS] <- as.numeric(phy$node.label)
	}
	
	model.set.final<-rate.mat.set.paint(phy,data.sort,nregimes,ntraits,model)
	
	if(!is.null(rate.mat)){
		rate <- rate.mat
		rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
		model.set.final$rate <- rate
		model.set.final$index.matrix <- rate.mat
	}

	lower = rep(lb, model.set.final$np)
	upper = rep(ub, model.set.final$np)
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	
	if(!is.null(p)){
		cat("Calculating likelihood from a set of fixed parameters", "\n")
		out<-NULL
		out$solution<-p
		out$objective<-dev.corpaint(out$solution,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
		loglik <- -out$objective
		est.pars<-out$solution
	}
	else{
		if(is.null(ip)){
			cat("Initializing...", "\n")
			#If the analysis is to be run a single processor:
			#Sets parameter settings for random restarts by taking the parsimony score and dividing
			#by the total length of the tree
			model.set.init<-rate.mat.set.paint(phy,data.sort,nregimes,ntraits,model="ER")
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
			dat<-as.matrix(data.sort)
			dat<-phyDat(dat,type="USER", levels=c("0","1"))
			par.score<-parsimony(phy, dat, method="fitch")
			tl <- sum(phy$edge.length)
			mean = par.score/tl
			ip<-rexp(1, 1/mean)
			lower.init = rep(lb, model.set.init$np)
			upper.init = rep(ub, model.set.init$np)
			init = nloptr(x0=rep(ip, length.out = model.set.init$np), eval_f=dev.corpaint, lb=lower.init, ub=upper.init, opts=opts, phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,regimes=regimes,root.p=root.p)
			cat("Finished. Beginning thorough search...", "\n")
			lower = rep(lb, model.set.final$np)
			upper = rep(ub, model.set.final$np)	
			out = nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.corpaint, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,regimes=regimes,root.p=root.p)
			loglik <- -out$objective
			est.pars<-out$solution
		}
		#If a user-specified starting value(s) is supplied:
		else{
			cat("Beginning subplex optimization routine -- Starting value(s):", ip, "\n")
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
			out = nloptr(x0=rep(ip, length.out = model.set.final$np), eval_f=dev.corpaint, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,regimes=regimes,root.p=root.p)
			loglik <- -out$objective
			est.pars<-out$solution
		}
	}
	#Starts the summarization process:
	cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
	
	lik.anc <- ancRECON.paint(phy, data, est.pars, hrm=FALSE, rate.cat=NULL, ntraits=ntraits, method=node.states, model=model, root.p=root.p)
	if(node.states == "marginal" || node.states == "scaled"){
		pr<-apply(lik.anc$lik.anc.states,1,which.max)
		phy$node.label <- pr
		tip.states <- NULL
	}
	if(ntraits==1){
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				colnames(lik.anc$lik.anc.states) <-  c("P(0)","P(1)")
			}
		}
	}
	if(ntraits==2){
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				colnames(lik.anc$lik.anc.states) <-  c("P(0,0)","P(0,1)","P(1,0)","P(1,1)")
			}
		}
	}
	if(ntraits==3){
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				colnames(lik.anc$lik.anc.states) <-  c("P(0,0,0)","P(1,0,0)","P(0,1,0)","P(0,0,1)","P(1,1,0)","P(1,0,1)","P(0,1,1)","P(1,1,1)")
			}
		}
	}
	if(node.states == "joint"){
		phy$node.label <- lik.anc$lik.anc.states
		tip.states <- lik.anc$lik.tip.states
	}

	cat("Finished. Performing diagnostic tests.", "\n")
	
	#Approximates the Hessian using the numDeriv function
	if(diagn==TRUE){
		h <- hessian(func=dev.corpaint, x=est.pars, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,regimes=regimes,root.p=root.p)
		hess.eig <- eigen(h,symmetric=TRUE)
		eigval<-signif(hess.eig$values,2)
		eigvect<-round(hess.eig$vectors, 2)
	}
	else{
		h=matrix(0,dim(model.set.final$index.matrix[[1]])[1],dim(model.set.final$index.matrix[[1]])[1])
		eigval=NULL
		eigvect=NULL
	}
	solution<-solution.se<-model.set.final$index.matrix
	for(i in 1:nregimes){
		solution[[i]] <- matrix(est.pars[model.set.final$index.matrix[[i]]], dim(model.set.final$index.matrix[[i]]))
		solution.se[[i]] <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix[[i]]], dim(model.set.final$index.matrix[[i]]))
		if(ntraits==1){
			rownames(solution[[i]]) <- rownames(solution.se[[i]]) <- c("(0)","(1)")
			colnames(solution[[i]]) <- colnames(solution.se[[i]]) <- c("(0)","(1)")
		}
		if(ntraits==2){
			rownames(solution[[i]]) <- rownames(solution.se[[i]]) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
			colnames(solution[[i]]) <- colnames(solution.se[[i]]) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
		}
		if(ntraits==3){
			rownames(solution[[i]]) <- rownames(solution.se[[i]]) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
			colnames(solution[[i]]) <- colnames(solution.se[[i]]) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
		}
	}

	obj = list(loglik = loglik, AIC = -2*loglik+2*model.set.final$np,AICc = -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))),ntraits=ntraits, solution=solution, solution.se=solution.se, index.mat=model.set.final$index.matrix, opts=opts, data=data.sort, phy=phy, states=lik.anc$lik.anc.states, tip.states=tip.states, iterations=out$iterations, eigval=eigval, eigvect=eigvect) 
	class(obj)<-"corpaint"
	return(obj)
}

#Print function
print.corpaint<-function(x,...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$ntraits,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","N.traits","ntax")
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

dev.corpaint<-function(p,phy,liks,Q,rate,regimes,root.p){
	
	nregimes<-length(levels(factor(regimes)))
	root.reg<-as.numeric(phy$node.label[1])
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	TIPS <- 1:nb.tip
	comp <- numeric(nb.tip + nb.node)
	phy <- reorder(phy, "pruningwise")
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])

	if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)
	
	for(i in 1:nregimes){
		Q[[i]]<- matrix(c(p, 0)[rate[[i]]],dim(rate[[1]])[1],dim(rate[[1]])[2])
		diag(Q[[i]]) <- -rowSums(Q[[i]])	
	}
	for (i  in seq(from = 1, length.out = nb.node)){
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		v <- 1
		for (desIndex in sequence(length(desRows))){
			v<-v*expm(Q[[regimes[focal]]] * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
		}
		comp[focal] <- sum(v)
		liks[focal, ] <- v/comp[focal]
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

rate.mat.set.paint<-function(phy,data.sort,nregimes,ntraits,model){

	k=ntraits
	nl=2
	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	if(ntraits==1){
		factored <- factorData(data.sort) # was acting on data, not data.sort			
		nl <- ncol(factored)
		obj <- NULL
		nb.tip<-length(phy$tip.label)
		nb.node <- phy$Nnode
		if (is.character(model)) {
			
			if(model=="ER"){
				rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=1, nstates=2, model="ER")
				np <- nregimes
				for(i in 2:nregimes){
					rate[[i]]<-rate[[i]]+i-1
				}
				index.matrix<-rate
				for(i in 1:nregimes){
					rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
				}
			}
			
			if(model=="ARD"){
				rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=1, model="ARD")
				np <- length(rate[[1]][is.na(rate[[1]])==TRUE])*nregimes
				
				for(i in 2:nregimes){
					rate[[i]]<-rate[[i]]+(length(rate[[1]][is.na(rate[[1]])==TRUE])*(i-1))
				}
				index.matrix<-rate
				for(i in 1:nregimes){
					rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
				}
			}
		}
		
		stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
		for(column in 1:nl){
			stateTable <- cbind(stateTable,factored[,column])
		}
		colnames(stateTable) <- colnames(factored)
		
		ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
		liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
		rownames(liks) <- NULL
	}	
	if(ntraits==2){
		if (is.character(model)) {

			if(model=="ER"){
				rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=2, model="ER")
				np <- nregimes
				for(i in 2:nregimes){
					rate[[i]]<-rate[[i]]+i-1
				}
				index.matrix<-rate
				for(i in 1:nregimes){
					rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
				}
			}
			if(model=="ARD"){
				rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=2, model="ARD")
				np <- length(rate[[1]][is.na(rate[[1]])==TRUE])*nregimes
				
				for(i in 2:nregimes){
					rate[[i]]<-rate[[i]]+(length(rate[[1]][is.na(rate[[1]])==TRUE])*(i-1))
				}
				index.matrix<-rate
				for(i in 1:nregimes){
					rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
				}
			}
		}	
		
		x<-data.sort[,1]
		y<-data.sort[,2]
		
		liks <- matrix(0, nb.tip + nb.node, nl^k)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(is.na(x[i])){
				x[i]=2
				y[i]=2
			}
		}
		for(i in 1:nb.tip){
			if(x[i]==0 & y[i]==0){liks[i,1]=1}
			if(x[i]==0 & y[i]==1){liks[i,2]=1}
			if(x[i]==1 & y[i]==0){liks[i,3]=1}
			if(x[i]==1 & y[i]==1){liks[i,4]=1}
			if(x[i]==2 & y[i]==2){liks[i,1:4]=1}
		}
	}
	if(ntraits==3){
		rate<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
		index.matrix<-rate
		rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
		
		x<-data.sort[,1]
		y<-data.sort[,2]
		z<-data.sort[,3]
		
		liks <- matrix(0, nb.tip + nb.node, nl^k)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(is.na(x[i])){
				x[i]=2
				y[i]=2
				z[i]=2
			}
		}
		for(i in 1:nb.tip){
			if(x[i]==0 & y[i]==0 & z[i]==0){liks[i,1]=1}
			if(x[i]==1 & y[i]==0 & z[i]==0){liks[i,2]=1}
			if(x[i]==0 & y[i]==1 & z[i]==0){liks[i,3]=1}
			if(x[i]==0 & y[i]==0 & z[i]==1){liks[i,4]=1}
			if(x[i]==1 & y[i]==1 & z[i]==0){liks[i,5]=1}
			if(x[i]==1 & y[i]==0 & z[i]==1){liks[i,6]=1}
			if(x[i]==0 & y[i]==1 & z[i]==1){liks[i,7]=1}
			if(x[i]==1 & y[i]==1 & z[i]==1){liks[i,8]=1}
			if(x[i]==2 & y[i]==2 & z[i]==2){liks[i,1:8]=1}
		}
	}
	Q <- rate
	for(i in 1:nregimes){
		Q[[i]]<-Q[[i]]*0
	}
	obj$np<-np
	obj$rate<-rate
	obj$index.matrix<-index.matrix
	obj$liks<-liks
	obj$Q<-Q
	
	obj
}

ancRECON.paint <- function(phy, data, p, method=c("joint", "marginal", "scaled"), hrm=TRUE, rate.cat, ntraits=NULL, charnum=NULL, rate.mat=NULL, model=c("ER", "SYM", "ARD"), root.p=NULL){
	
	#Note: Does not like zero branches at the tips. Here I extend these branches by just a bit:
	phy$edge.length[phy$edge.length<=1e-5]=1e-5

	if(hrm==FALSE){
		if(ntraits==1){
			data.sort<-data.frame(data[,2],data[,3],row.names=data[,1])
		}
		if(ntraits==2){
			data.sort<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
	}
		if(ntraits==3){
			data.sort<-data.frame(data[,2], data[,3], data[,4], data[,5], row.names=data[,1])
		}
	}
	else{
		data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
	}
	data.sort<-data.sort[phy$tip.label,]
	root.reg<-as.numeric(phy$node.label[1])

	#Some initial values for use later
	k=2
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	
	root.p=root.p

	if(ntraits==1){
		tot.states <- factor(c(phy$node.label,as.character(data[,3])))
		nregimes <- length(levels(tot.states))
		regimes <- numeric(nb.tip+nb.node)
		TIPS <- 1:nb.tip		
		regimes[TIPS]<-data[,3]
		regimes[-TIPS]<-as.numeric(phy$node.label)
	}
	if(ntraits==2){
		tot.states <- factor(c(phy$node.label,as.character(data[,4])))
		nregimes <- length(levels(tot.states))
		regimes <- numeric(nb.tip+nb.node)
		TIPS <- 1:nb.tip		
		regimes[TIPS]<-data[,4]
		regimes[-TIPS]<-as.numeric(phy$node.label)
	}
	if(ntraits==3){
		tot.states <- factor(c(phy$node.label,as.character(data[,5])))
		nregimes <- length(levels(tot.states))
		regimes <- numeric(nb.tip+nb.node)
		TIPS <- 1:nb.tip				
		regimes[TIPS] <- data[,5]
		regimes[-TIPS] <- as.numeric(phy$node.label)
	}

	if(hrm==FALSE){
		#Imported from Jeffs rayDISC -- will clean up later, but for now, it works fine:
		if(ntraits==1){
			k <- 1
			factored <- factorData(data.sort) # was acting on data, not data.sort			
			nl <- ncol(factored)
			nb.tip<-length(phy$tip.label)
			nb.node <- phy$Nnode
			if (is.character(model)) {
				
				if(model=="ER"){
					rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=1, nstates=nl, model="ER")
					np <- nregimes
					for(i in 2:nregimes){
						rate[[i]]<-rate[[i]]+i-1
					}
					index.matrix<-rate
					for(i in 1:nregimes){
						rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
					}
				}
				if(model=="ARD"){
					rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=1, nstates=nl, model="ARD")
					np <- length(rate[[1]][is.na(rate[[1]])==TRUE])*nregimes
					
					for(i in 2:nregimes){
						rate[[i]]<-rate[[i]]+(length(rate[[1]][is.na(rate[[1]])==TRUE])*(i-1))
					}
					index.matrix<-rate
					for(i in 1:nregimes){
						rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
					}
				}
			}

			stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
			for(column in 1:nl){
				stateTable <- cbind(stateTable,factored[,column])
			}
			colnames(stateTable) <- colnames(factored)
			
			ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
			liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
			rownames(liks) <- NULL
		}
		if(ntraits==2){
			k=2
			nl=2
			if(is.null(rate.mat)){
				if(model=="ER"){
					rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=2, model="ER")
					np <- nregimes
					for(i in 2:nregimes){
						rate[[i]]<-rate[[i]]+i-1
					}
					index.matrix<-rate
					for(i in 1:nregimes){
						rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
					}
				}
				if(model=="ARD"){
					rate <- lapply(1:nregimes,rate.mat.maker, hrm=FALSE, ntraits=2, model="ARD")
					np <- length(rate[[1]][is.na(rate[[1]])==TRUE])*nregimes
					
					for(i in 2:nregimes){
						rate[[i]]<-rate[[i]]+(length(rate[[1]][is.na(rate[[1]])==TRUE])*(i-1))
					}
					index.matrix<-rate
					for(i in 1:nregimes){
						rate[[i]][is.na(rate[[i]])]<-max(rate[[nregimes]],na.rm=TRUE)+1
					}
				}
			}
			else{
				rate<-rate.mat
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			x<-data.sort[,1]
			y<-data.sort[,2]
			
			liks <- matrix(0, nb.tip + nb.node, nl^k)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(is.na(x[i])){
					x[i]=2
				}
				if(is.na(y[i])){
					y[i]=2
				}
			}
			for(i in 1:nb.tip){
				if(x[i]==0 & y[i]==0){liks[i,1]=1}
				if(x[i]==0 & y[i]==1){liks[i,2]=1}
				if(x[i]==1 & y[i]==0){liks[i,3]=1}
				if(x[i]==1 & y[i]==1){liks[i,4]=1}
				if(x[i]==2 & y[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==2 & y[i]==1){liks[i,c(2,4)]=1}
				if(x[i]==0 & y[i]==2){liks[i,c(1,2)]=1}
				if(x[i]==1 & y[i]==2){liks[i,c(3,4)]=1}
				if(x[i]==2 & y[i]==2){liks[i,1:4]=1}
			}
		}
		if(ntraits==3){
			k=3
			nl=2
			if(is.null(rate.mat)){
				rate<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			else{
				rate<-rate.mat
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}			
			x<-data.sort[,1]
			y<-data.sort[,2]
			z<-data.sort[,3]
			
			liks <- matrix(0, nb.tip + nb.node, nl^k)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(is.na(x[i])){
					x[i]=2
				}
				if(is.na(y[i])){
					y[i]=2
				}
				if(is.na(z[i])){
					z[i]=2
				}
			}
			for(i in 1:nb.tip){
				if(x[i]==0 & y[i]==0 & z[i]==0){liks[i,1]=1}
				if(x[i]==1 & y[i]==0 & z[i]==0){liks[i,2]=1}
				if(x[i]==0 & y[i]==1 & z[i]==0){liks[i,3]=1}
				if(x[i]==0 & y[i]==0 & z[i]==1){liks[i,4]=1}
				if(x[i]==1 & y[i]==1 & z[i]==0){liks[i,5]=1}
				if(x[i]==1 & y[i]==0 & z[i]==1){liks[i,6]=1}
				if(x[i]==0 & y[i]==1 & z[i]==1){liks[i,7]=1}
				if(x[i]==1 & y[i]==1 & z[i]==1){liks[i,8]=1}
				#If x is ambiguous but the rest are not:
				if(x[i]==2 & y[i]==0 & z[i]==0){liks[i,c(1,2)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==0){liks[i,c(3,5)]=1}
				if(x[i]==2 & y[i]==0 & z[i]==1){liks[i,c(4,6)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==1){liks[i,c(7,8)]=1}
				#If y is ambiguous but the rest are not:
				if(x[i]==0 & y[i]==2 & z[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==0){liks[i,c(2,5)]=1}
				if(x[i]==0 & y[i]==2 & z[i]==1){liks[i,c(4,7)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==1){liks[i,c(6,8)]=1}
				#If z is ambiguous but the rest are not:
				if(x[i]==0 & y[i]==0 & z[i]==2){liks[i,c(1,4)]=1}
				if(x[i]==0 & y[i]==1 & z[i]==2){liks[i,c(3,7)]=1}
				if(x[i]==1 & y[i]==0 & z[i]==2){liks[i,c(2,6)]=1}
				if(x[i]==1 & y[i]==1 & z[i]==2){liks[i,c(5,8)]=1}
				#If x and y is ambiguous but z is not:
				if(x[i]==2 & y[i]==2 & z[i]==0){liks[i,c(1,2,3,5)]=1}
				if(x[i]==2 & y[i]==2 & z[i]==1){liks[i,c(4,6,7,8)]=1}
				#If x and z is ambiguous but y is not:
				if(x[i]==2 & y[i]==0 & z[i]==2){liks[i,c(1,2,4,6)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==2){liks[i,c(3,5,7,8)]=1}
				#If y and z is ambiguous but x is not:
				if(x[i]==0 & y[i]==2 & z[i]==2){liks[i,c(1,3,4,7)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==2){liks[i,c(2,5,6,8)]=1}
				#All states are ambiguous:			
				if(x[i]==2 & y[i]==2 & z[i]==2){liks[i,1:8]=1}
			}
		}
		Q <- rate
		for(i in 1:nregimes){
			Q[[i]]<-Q[[i]]*0
		}
	}
	for(i in 1:nregimes){
		Q[[i]]<- matrix(c(p, 0)[rate[[i]]],dim(rate[[1]])[1],dim(rate[[1]])[2])
		diag(Q[[i]]) <- -rowSums(Q[[i]])	
	}
	phy <- reorder(phy, "pruningwise")
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])

	if(method=="joint"){
		lik.states<-numeric(nb.tip + nb.node)
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		for (i  in seq(from = 1, length.out = nb.node)) {
			#The ancestral node at row i is called focal:
			focal <- anc[i]
			#Get descendant information of focal:
			desRows<-which(phy$edge[,1]==focal)
			#Get node information for each descendant:
			desNodes<-phy$edge[desRows,2]
			#Initiates a loop to check if any nodes are tips:
			for (desIndex in sequence(length(desRows))){
				#If a tip calculate C_y(i) for the tips and stores in liks matrix:
				if(any(desNodes[desIndex]==phy$edge[,1])==FALSE){
					liks[desNodes[desIndex],] <- expm(Q[[regimes[focal]]] * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
					#Divide by the sum of the liks to deal with underflow issues:
					liks[desNodes[desIndex],] <- liks[desNodes[desIndex],]/sum(liks[desNodes[desIndex],])
					#Collects the likeliest state at the tips:
					comp[desNodes[desIndex],] <- which.max(liks[desNodes[desIndex],])
				}
			}
			#Collects t_z, or the branch subtending focal:
			tz<-phy$edge.length[which(phy$edge[,2] == focal)]	
			if(length(tz)==0){
				#The focal node is the root, calculate P_k:
				if(is.null(root.p)){
					root.state=1
					for (desIndex in sequence(length(desRows))){
						#This is the basic marginal calculation:
						root.state <- root.state * expm(Q[[regimes[focal]]] * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
					}
					if(is.null(root.p)){
						liks[focal, ] <- root.state
					}
					else{
						if(is.character(root.p)){
							equil.root <- NULL
							for(i in 1:ncol(Q[[root.reg]])){
								posrows <- which(Q[[root.reg]][,i] >= 0)
								rowsum <- sum(Q[[root.reg]][posrows,i])
								poscols <- which(Q[[root.reg]][i,] >= 0)
								colsum <- sum(Q[[root.reg]][i,poscols])
								equil.root <- c(equil.root,rowsum/(rowsum+colsum))
							}
							liks[focal, ] <- root.state * equil.root
						}
						else{
							liks[focal, ] <- root.p
						}
					}
					liks[focal, ] <- liks[focal,] / sum(liks[focal,])
				}
			}
			else{
				#Calculates P_ij(t_z):
				Pij <- expm(Q[[regimes[focal]]] * tz, method=c("Ward77"))
				#Calculates L_z(i):
				if(hrm==TRUE){
					v<-c(rep(1, k*rate.cat))
				}
				if(hrm==FALSE){
					v<-c(rep(1, nl^k))
				}
				for (desIndex in sequence(length(desRows))){
					v = v * liks[desNodes[desIndex],]
				}
				#Finishes L_z(i):
				L <- t(Pij) * v
				#Collects which is the highest likelihood and which state it corresponds to:
				liks[focal,] <- apply(L, 2, max)
				comp[focal,] <- apply(L, 2, which.max)
				#Divide by the sum of the liks to deal with underflow issues:
				liks[focal,] <- liks[focal,]/sum(liks[focal,])
			}
		}
		root <- nb.tip + 1L
		lik.states[root] <- which.max(liks[root,])
		N <- dim(phy$edge)[1]
		for(i in N:1){
			des <- phy$edge[i,2]
			tmp <- which.max(liks[des,])
			lik.states[des] <- comp[des,tmp]
		}
		#Outputs likeliest tip states
		obj$lik.tip.states <- lik.states[TIPS]
		#Outputs likeliest node states
		obj$lik.anc.states <- lik.states[-TIPS]
	}
	if(method=="marginal"){
		#A temporary likelihood matrix so that the original does not get written over:
		liks.down<-liks
		#root equilibrium frequencies
		if(is.null(root.p)){
			equil.root<-rep(1/dim(Q[[1]])[2], dim(Q[[1]])[2])
		}
		else{			
			if(is.character(root.p)){
				equil.root <- NULL
				for(i in 1:ncol(Q[[root.reg]])){
					posrows <- which(Q[[root.reg]][,i] >= 0)
					rowsum <- sum(Q[[root.reg]][posrows,i])
					poscols <- which(Q[[root.reg]][i,] >= 0)
					colsum <- sum(Q[[root.reg]][i,poscols])
					equil.root <- c(equil.root,rowsum/(rowsum+colsum))
				}
			}
			else{
				equil.root=root.p
			}
		}
		#A transpose of Q for assessing probability of j to i, rather than i to j:
		tranQ<-Q		
		for(i in 1:nregimes){		
			tranQ[[i]]<-t(Q[[i]])
		}
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		#The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
		for (i  in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q[[regimes[focal]]] * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
			}
			comp[focal] <- sum(v)
			liks.down[focal, ] <- v/comp[focal]
		}
		root <- nb.tip + 1L
		#Enter the root defined root probabilities if they are supplied by the user:
		if(!is.null(root.p)){
			root <- nb.tip + 1L	
			liks.down[root, ]<-root.p
		}
		#The up-pass 
		liks.up<-liks
		states<-apply(liks,1,which.max)
		N <- dim(phy$edge)[1]
		comp <- numeric(nb.tip + nb.node)
		for(i in length(anc):1){
			focal <- anc[i]
			if(!focal==root){
				#Gets mother and sister information of focal:			
				focalRow<-which(phy$edge[,2]==focal)
				motherRow<-which(phy$edge[,1]==phy$edge[focalRow,1])
				motherNode<-phy$edge[focalRow,1]
				desNodes<-phy$edge[motherRow,2]
				sisterNodes<-desNodes[(which(!desNodes==focal))]
				sisterRows<-which(phy$edge[,2]%in%sisterNodes==TRUE)		
				#If the mother is not the root then you are calculating the probability of the being in either state.
				#But note we are assessing the reverse transition, j to i, rather than i to j, so we transpose Q to carry out this calculation:
				if(motherNode!=root){
					v <- expm(tranQ[[regimes[motherNode]]] * phy$edge.length[which(phy$edge[,2]==motherNode)], method=c("Ward77")) %*% liks.up[motherNode,]		
				}
				#If the mother is the root then just use the marginal. This can also be the prior, which I think is the equilibrium frequency. 
				#But for now we are just going to use the marginal at the root -- it is unclear what Mesquite does.
				else{
					v <- equil.root
				}
				#Now calculate the probability that each sister is in either state. Sister can be more than 1 when the node is a polytomy. 
				#This is essentially calculating the product of the mothers probability and the sisters probability:
				for (sisterIndex in sequence(length(sisterRows))){
					v <- v*expm(Q[[regimes[sisterNodes[sisterIndex]]]] * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
				}
				comp[focal] <- sum(v)
				liks.up[focal,] <- v/comp[focal]
			}
		}

		#The final pass
		liks.final<-liks
		comp <- numeric(nb.tip + nb.node)
		for (i  in seq(from = 1, length.out = nb.node-1)) { # In this final pass, root is never encountered.  But its OK, because root likelihoods are set after loop.
			#the ancestral node at row i is called focal
			focal <- anc[i]
			focalRow<-which(phy$edge[,2]==focal)
			motherRow<-which(phy$edge[,1]==phy$edge[focalRow,1])
			motherNode<-phy$edge[focalRow,1]
			#Now you are assessing the change along the branch subtending the focal by multiplying the probability of 
			#everything at and above focal by the probability of the mother and all the sisters given time t:
			v <- liks.down[focal,]*expm(tranQ[[regimes[motherNode]]] * phy$edge.length[focalRow], method=c("Ward77")) %*% liks.up[focal,]
			comp[focal] <- sum(v)
			liks.final[focal, ] <- v/comp[focal]
		}
		#Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
		liks.final[root,] <- liks.down[root,] * equil.root
		
		root.final <- liks.down[root,] * equil.root
		comproot <- sum(root.final)
		liks.final[root,] <- root.final/comproot
		#Reports just the probabilities at internal nodes:
		obj$lik.anc.states <- liks.final[-TIPS, ]
	}	
	
	if(method=="scaled"){
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		#The same algorithm as in the main function. See comments in either corHMM.R or corDISC.R for details:
		for (i  in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q[[regimes[focal]]] * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
			}
			comp[focal] <- sum(v)
			liks[focal, ] <- v/comp[focal]
		}
		if(!is.null(root.p)){
			root <- nb.tip + 1L	
			if(is.character(root.p)){
				equil.root <- NULL
				for(i in 1:ncol(Q[[root.reg]])){
					posrows <- which(Q[[root.reg]][,i] >= 0)
					rowsum <- sum(Q[[root.reg]][posrows,i])
					poscols <- which(Q[[root.reg]][i,] >= 0)
					colsum <- sum(Q[[root.reg]][i,poscols])
					equil.root <- c(equil.root,rowsum/(rowsum+colsum))
				}
				liks[root, ] <- (liks[root,] * equil.root) / sum((liks[root,] * equil.root))
			}
			else{
				liks[root, ] <- root.p
			}
		}
		obj$lik.anc.states <- liks[-TIPS, ]
	}	
	obj
}


