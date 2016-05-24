#OUwie Master Controller

#written by Jeremy M. Beaulieu

#Fits the Ornstein-Uhlenbeck model of continuous characters evolving under discrete selective 
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels 
#and a trait file. The trait file must be in the following order: Species names, Regime, and 
#continuous trait. Different models can be specified -- Brownian motion (BM), multiple rate BM (BMS)
#global OU (OU1), multiple regime OU (OUM), multiple sigmas (OUMV), multiple alphas (OUMA), 
#and the multiple alphas and sigmas (OUMVA). 

OUwie<-function(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA", "TrendyM", "TrendyMS"), simmap.tree=FALSE, scaleHeight=FALSE, root.station=TRUE, clade=NULL, mserr="none", diagn=FALSE, quiet=FALSE, warn=TRUE){
	#Makes sure the data is in the same order as the tip labels
	if(mserr=="none" | mserr=="est"){
		data<-data.frame(data[,2], data[,3], row.names=data[,1])
		data<-data[phy$tip.label,]
	}
	if(mserr=="known"){
		if(!dim(data)[2]==4){
			stop("You specified measurement error should be incorporated, but this information is missing")
		}
		else{
			data<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
			data<-data[phy$tip.label,]
		}
	}
	#Values to be used throughout
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	#Will label the clade of interest if the user so chooses:
	if(is.null(clade)){
		phy=phy
	}
	if(!is.null(clade) & simmap.tree==FALSE){
		node<-mrca(phy)[clade[1],clade[2]]
		int<-c(node,Descendants(phy,node, "all"))
		tips<-int[int<ntips]
		data[,1]<-1
		data[tips,1]<-2
		int<-int[int>ntips]
		phy$node.label<-rep(1,phy$Nnode)
		pp<-int-length(phy$tip.label)
		phy$node.label[pp]<-2
	}
	if (is.character(model)) {
		if (model == "BM1"| model == "OU1"){
			simmap.tree=FALSE
		}
	}
	if(simmap.tree==TRUE){
		k<-length(colnames(phy$mapped.edge))
		tot.states<-factor(colnames(phy$mapped.edge))
		tip.states<-factor(data[,1])
		data[,1]<-as.numeric(tip.states)
		
		#A boolean for whether the root theta should be estimated -- default is that it should be.
		root.station=root.station
		#Obtains the state at the root
		root.state=which(colnames(phy$mapped.edge)==names(phy$maps[[1]][1]))
		##Begins the construction of the edges matrix -- similar to the ouch format##
		edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
		
		if(scaleHeight==TRUE){
			edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
		}
		edges=edges[sort.list(edges[,3]),]
		
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
	}
	if(simmap.tree==FALSE){
		#Obtain a a list of all the regime states. This is a solution for instances when tip states and 
		#the internal nodes are not of equal length:
		tot.states<-factor(c(phy$node.label,as.character(data[,1])))
		k<-length(levels(tot.states))
		int.states<-factor(phy$node.label)
		phy$node.label=as.numeric(int.states)		
		tip.states<-factor(data[,1])
		data[,1]<-as.numeric(tip.states)
		
		#A boolean for whether the root theta should be estimated -- default is that it should be.
		root.station=root.station
		if (is.character(model)) {			
			if (model == "BM1"| model == "OU1"){
				##Begins the construction of the edges matrix -- similar to the ouch format##
				#Makes a vector of absolute times in proportion of the total length of the tree
				phy$node.label<-sample(c(1:k), phy$Nnode, replace=T)
				int.states=length(levels(tip.states))
				#Since we only really have one global regime, make up the internal nodes -- this could be improved
				phy$node.label<-as.numeric(int.states)				
				#Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
				root.state<-1
				#New tree matrix to be used for subsetting regimes
				edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
				if(scaleHeight==TRUE){
					edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
				}
				edges=edges[sort.list(edges[,3]),]
				
				regime <- matrix(0,nrow=length(edges[,1]),ncol=k)
				regime[,1]<-1
				if (k>1) {
					regime[,2:k]<-0
				}
				edges=cbind(edges,regime)
			}
			else{			
				#Obtain root state and internal node labels
				root.state<-phy$node.label[1]
				int.state<-phy$node.label[-1]
				#New tree matrix to be used for subsetting regimes
				edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
				if(scaleHeight==TRUE){
					edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
				}
				edges=edges[sort.list(edges[,3]),]
				
				mm<-c(data[,1],int.state)
				regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
				#Generates an indicator matrix from the regime vector
				for (i in 1:length(mm)) {
					regime[i,mm[i]] <- 1 
				}
				#Finishes the edges matrix
				edges=cbind(edges,regime)
			}
		}
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
	}
	
	#Initializes an integer specifying whether or not the trendy model is to be invoked. Assume 0, unless otherwise stated later:
	trendy=0
	#Data:
	x<-as.matrix(data[,2])
	
	#Matches the model with the appropriate parameter matrix structure
	if (is.character(model)) {
		index.mat<-matrix(0,2,k)
		Rate.mat <- matrix(1, 2, k)
		if (model == "BM1"){
			np=1
			index<-matrix(TRUE,2,k)
			if(mserr=="est"){
				index.mat[1,1:k]<-np+2
			}
			else{
				index.mat[1,1:k]<-np+1				
			}
			index.mat[2,1:k]<-1
			param.count<-np+1
			bool=TRUE
		}
		#The group mean model of Thomas et al, is trivial: set bool to be TRUE:
		if (model == "BMS"){
			np=k
			index<-matrix(TRUE,2,k)
			if(mserr=="est"){
				index.mat[1,1:k]<-np+2
			}
			else{
				index.mat[1,1:k]<-np+1
			}
			index.mat[2,1:k]<-1:np
#			if(root.station==TRUE){
#				param.count<-np+k
#			}
#			if(root.station==FALSE){
			param.count<-np+1
#			}			
			bool=FALSE
		}
		if (model == "OU1"){
			np=2
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1
			index.mat[2,1:k]<-2
			if(root.station==TRUE){
				param.count<-np+1
			}
			if(root.station==FALSE){
				param.count<-np+2
			}
			bool=root.station
		}
		if (model == "OUM"){
			np=2
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1
			index.mat[2,1:k]<-2
			if(root.station==TRUE){
				param.count<-np+k
			}
			if(root.station==FALSE){
				param.count<-np+k+1
			}
			bool=root.station
		}
		if (model == "OUMV") {
			np=k+1
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1
			index.mat[2,1:k]<-2:(k+1)
			if(root.station==TRUE){
				param.count<-np+k
			}
			if(root.station==FALSE){
				param.count<-np+k+1
			}
			bool=root.station
		}
		if (model == "OUMA") {
			np=k+1
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1:k
			index.mat[2,1:k]<-k+1
			if(root.station==TRUE){
				param.count<-np+k
			}
			if(root.station==FALSE){
				param.count<-np+k+1
			}
			bool=root.station
		}
		if (model == "OUMVA") {
			np=k*2
			index<-matrix(TRUE,2,k)
			index.mat[index]<-1:(k*2)
			if(root.station==TRUE){
				param.count<-np+k
			}
			if(root.station==FALSE){
				param.count<-np+k+1
			}
			bool=root.station
		}
		if (model == "TrendyM"){
			index.mat<-matrix(0,3,k)
			Rate.mat <- matrix(1, 3, k)
			np=k+2 #We have k regimes, but 1 sigma, plus the root value
			index<-matrix(TRUE,3,k)
			if(mserr=="est"){
				index.mat[1,1:k]<-np+2
			}
			else{
				index.mat[1,1:k]<-np+1				
			}
			index.mat[2,1:k]<-1
			index.mat[3,1:k]<-2:(np-1)
			param.count<-np
			bool=TRUE
			#trendy=1
		}
		if (model == "TrendyMS"){
			index.mat<-matrix(0,3,k)
			Rate.mat <- matrix(1, 3, k)
			np= (k*2)+1 #We have k regimes and assume k trends and k sigmas, plus the root value
			index<-matrix(TRUE,3,k)
			if(mserr=="est"){
				index.mat[1,1:k]<-np+2
			}
			else{
				index.mat[1,1:k]<-np+1
			}
			index.mat[2,1:k]<-1:k
			index.mat[3,1:k]<-(k+1):(np-1)
			param.count<-np
			bool=FALSE
			#trendy=1
		}
	}

	if(warn==TRUE){
		if(param.count > (ntips/10)){
			warning("You might not have enough data to fit this model well", call.=FALSE, immediate.=TRUE)
		}
	}
	#Likelihood function for estimating model parameters
	dev<-function(p, index.mat, edges, mserr, trendy){
				
		p <- exp(p)
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		N<-length(x[,1])
		root.par.index=length(p)
		V<-varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight)
		if (any(is.nan(diag(V))) || any(is.infinite(diag(V)))) return(1000000)		
		if(mserr=="known"){
			diag(V)<-diag(V)+(data[,3]^2)
		}
		if(mserr=="est"){
			root.par.index=length(p) - 1
			diag(V)<-diag(V)+(p[length(p)])
		}
		if(trendy == 0){
			W<-weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, assume.station=bool)
			theta<-Inf
			try(theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x, silent=TRUE)
			if(any(theta==Inf)){
				return(10000000)
			}
			DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
			#However, sometimes this fails (not sure yet why) so I just toggle between this and another approach:
			if(!is.finite(DET)){
				DET <- determinant(V, logarithm=TRUE)
				logl <- -.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
			}else{
				logl <- -.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET)-.5*(N*log(2*pi))
			}
		}else{
			E_a <- Inf
			try(E_a <- expected.trendy(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, root.value=p[root.par.index]))
			if(any(E_a==Inf)){
				return(10000000)
			}
			DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
			#However, sometimes this fails (not sure yet why) so I just toggle between this and another approach:
			if(!is.finite(DET)){
				DET <- determinant(V, logarithm=TRUE)
				logl <- -.5*(t(E_a-x)%*%pseudoinverse(V)%*%(E_a-x))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
			}else{
				logl <- -.5*(t(E_a-x)%*%pseudoinverse(V)%*%(E_a-x))-.5*as.numeric(DET)-.5*(N*log(2*pi))
			}
		}		
		#When the model includes alpha, the values of V can get too small, the modulus does not seem correct and the loglik becomes unstable. This is one solution:
		if(!is.finite(logl)){
			return(10000000)
		}
		return(-logl)
	}
	
	if(quiet==FALSE){
		cat("Initializing...", "\n")
	}

	lb = -20
	ub = 20
	lower = rep(-20, np)
	upper = rep(20, np)
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	
	if(model == "OU1" | model == "OUM" | model == "OUMV" | model == "OUMA" | model == "OUMVA"){
		if(scaleHeight==TRUE){
			d <- max(diag(vcv.phylo(phy)))
			phy$edge.length<-(phy$edge.length/d)
		}
		C.mat<-vcv.phylo(phy)
		a<-as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		sig<-as.numeric(t(x-a)%*%solve(C.mat)%*%(x-a)/n)
		init.np=2
		init.lower = rep(lb, init.np)
		init.upper = rep(ub, init.np)
		init.index.mat<-matrix(0,2,k)
		init.index.mat[1,1:k]<-1
		init.index.mat[2,1:k]<-2
		if(model == "OU1" | model == "OUM"){
			edges.tmp=edges
			if(dim(edges.tmp)[2] < 6){
				edges.tmp <- cbind(edges.tmp, rep(1, length(edges.tmp[,1])))
				edges.tmp <- cbind(edges.tmp,0)
			}else{
				edges.tmp[,6]<-rep(1, length(edges.tmp[,1]))
				edges.tmp[,7:(k+5)]<-0
			}
		}else{
			edges.tmp=edges
		}
		#Initial value for alpha is just the half life based on the entire length of the tree:
		init <- nloptr(x0=log(c(log(2)/max(branching.times(phy)),sig)), eval_f=dev, lb=init.lower, ub=init.upper, opts=opts, index.mat=init.index.mat, edges=edges.tmp, mserr=mserr, trendy=trendy)
		init.ip <- c(exp(init$solution[1]), exp(init$solution[2]))
		if(model=="OU1"){
			ip=init.ip
		}
		if(model=="OUMV" | model=="OUMA" | model=="OUM"){
			ip<-c(rep(init.ip[1],length(unique(index.mat[1,]))),rep(init.ip[2],length(unique(index.mat[2,]))))
		}
		if(model=="OUMVA"){
			ip<-c(rep(init.ip,k))
		}
		if(mserr=="est"){
			ip<-c(ip,sig)
			lower = c(lower,0)
			upper = c(upper,ub)
		}		
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		out = nloptr(x0=log(ip), eval_f=dev, lb=lower, ub=upper, opts=opts, index.mat=index.mat, edges=edges, mserr=mserr, trendy=trendy)
	}
	else{
		if(scaleHeight==TRUE){
			d <- max(diag(vcv.phylo(phy)))
			phy$edge.length<-(phy$edge.length/d)
		}
		#Starting values follow from phytools:
		C.mat<-vcv.phylo(phy)
		a<-as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		sig<-as.numeric(t(x-a)%*%solve(C.mat)%*%(x-a)/n)
		#####################
		if(model=="BMS"){
			ip=rep(sig,k)
		}
		else{
			if(model=="BM1"){
				ip=sig
			}
			if(model=="TrendyM"){
				#We assume that the starting trend values are zero:
				ip=c(sig, rep(0, param.count-2), a)
				lower = c(lb,rep(-1e6, param.count-2), -1e6)
				upper = c(ub,rep(1e6, param.count-2), 1e6)
			}
			if(model == "TrendyMS"){
				#We assume that the starting trend values are zero:
				ip=c(rep(sig, k), rep(0, (np-1) - k), a)
				lower = c(rep(lb,k),rep(-1e6, (np-1) - k), -1e6)
				upper = c(rep(ub,k),rep(1e6, (np-1) - k), 1e6)
			}
		}
		if(mserr=="est"){
			ip<-c(ip,sig)
			lower = c(lower,0)
			upper = c(upper,10)
		}
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		
		out = nloptr(x0=log(ip), eval_f=dev, lb=lower, ub=upper, opts=opts, index.mat=index.mat, edges=edges, mserr=mserr, trendy=trendy)
	}
	loglik <- -out$objective
	out$solution = exp(out$solution)

	#Takes estimated parameters from dev and calculates theta for each regime:
	dev.theta<-function(p, index.mat, edges=edges, mserr=mserr){
		tmp<-NULL
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		N<-length(x[,1])
		V<-varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight)
		W<-weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, assume.station=bool)
		if(mserr=="known"){
			diag(V)<-diag(V)+(data[,3]^2)
		}
		if(mserr=="est"){
			diag(V)<-diag(V)+(p[length(p)])
		}
		theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
		#Calculates the hat matrix:
		#H.mat<-W%*%pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)
		#Standard error of theta -- uses pseudoinverse to overcome singularity issues
		se<-sqrt(diag(pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)))
		tmp$res<-W%*%theta-x
		#Joins the vector of thetas with the vector of standard errors into a 2 column matrix for easy extraction at the summary stage
		tmp$theta.est<-cbind(theta,se)
		#Returns final GLS solution
		tmp
	}
	
	#Informs the user that the summarization has begun, output model-dependent and dependent on whether the root theta is to be estimated
	if(quiet==FALSE){
		cat("Finished. Summarizing results.", "\n")	
	}
	
	if(model == "TrendyM" | model == "TrendyMS"){
		theta<-NULL
		root.est <- out$solution[length(out$solution)]
		theta$theta.est <- matrix(root.est, k+1, 2)
		theta$theta.est[,2] <- 0 
		theta$res <- expected.trendy(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, root.value=root.est) - x
	}else{
		theta <- dev.theta(out$solution, index.mat, edges, mserr)
	}
	
	#Calculates the Hessian for use in calculating standard errors and whether the maximum likelihood solution was found
	if(diagn==TRUE){
		h <- hessian(x=log(out$solution), func=dev, index.mat=index.mat, edges=edges, mserr=mserr, trendy=trendy)
		#Using the corpcor package here to overcome possible NAs with calculating the SE
		solution<-matrix(out$solution[index.mat], dim(index.mat))
		solution.se<-matrix(sqrt(diag(pseudoinverse(h)))[index.mat], dim(index.mat))
		rownames(solution) <- rownames(solution.se) <- rownames(index.mat) <- c("alpha","sigma.sq")
		if(simmap.tree==FALSE){
			colnames(solution) <- colnames(solution.se) <- levels(tot.states)
		}
		if(simmap.tree==TRUE){
			colnames(solution) <- colnames(solution.se) <- c(colnames(phy$mapped.edge))
		}				
		#Eigendecomposition of the Hessian to assess reliability of likelihood estimates
		hess.eig<-eigen(h,symmetric=TRUE)
		#If eigenvect is TRUE then the eigenvector and index matrix will appear in the list of objects 
		eigval<-signif(hess.eig$values,2)
		eigvect<-round(hess.eig$vectors, 2)
		if(mserr=="est"){
			mserr.est<-out$solution[length(out$solution)]
			param.count<-param.count+1
		}
		else{
			mserr.est<-NULL
		}		
		obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))),model=model,solution=solution, mserr.est=mserr.est, theta=theta$theta.est, solution.se=solution.se, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, opts=opts, data=data, phy=phy, root.station=root.station, lb=lower, ub=upper, iterations=out$iterations, res=theta$res, eigval=eigval, eigvect=eigvect) 
		
	}
	if(diagn==FALSE){
		solution<-matrix(out$solution[index.mat], dim(index.mat))
		if(model=="TrendyM" | model=="TrendyMS"){
			rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq","trend")
		}else{
			rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq")
		}
		if(simmap.tree==FALSE){
			colnames(solution) <- levels(tot.states)
		}
		if(simmap.tree==TRUE){
			colnames(solution) <- c(colnames(phy$mapped.edge))
		}
		if(mserr=="est"){
			mserr.est<-out$solution[length(out$solution)]
			param.count<-param.count+1
		}
		else{
			mserr.est<-NULL
		}
		obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))),model=model,solution=solution, mserr.est=mserr.est, theta=theta$theta.est, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, opts=opts, data=data, phy=phy, root.station=root.station, lb=lower, ub=upper, iterations=out$iterations, res=theta$res) 
	}
	class(obj)<-"OUwie"		
	return(obj)
}

print.OUwie<-function(x, ...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$model,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","model","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	if (is.character(x$model)) {
		if (x$model == "BM1" | x$model == "BMS" | x$model == "TrendyM" | x$model == "TrendyMS"){
			param.est <- x$solution
#			if(x$root.station==FALSE){
			theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
#			}
#			else{
#				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
#			}
			rownames(theta.mat)<-c("estimate", "se")
			if(x$simmap.tree==FALSE){
				colnames(theta.mat) <- levels(x$tot.states)
			}
			if(x$simmap.tree==TRUE){
				colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
			}
			cat("Rates\n")
			print(param.est)
			cat("\n")
			cat("Optima\n")
			print(theta.mat)
			cat("\n")
		}
		if (x$root.station == TRUE | x$root.station==FALSE){
			if (x$model == "OU1"){
				param.est<- x$solution
				theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat) <- levels(x$tot.states)
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
				}
				cat("Rates\n")
				print(param.est)
				cat("\n")
				cat("\nOptima\n")
				print(theta.mat)
				cat("\n")
			}
		}
		if (x$root.station == TRUE){
			if (x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
				param.est<- x$solution
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat)<- levels(x$tot.states)
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
				}
				cat("\nRates\n")
				print(param.est)
				cat("\n")
				cat("Optima\n")
				print(theta.mat)
				cat("\n")
			}
		}
		if (x$root.station == FALSE){
			if (x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){ 
				param.est<- x$solution
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states))+1)
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat)<-c("root", levels(x$tot.states))
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat)<-c("root", colnames(x$phy$mapped.edge))
				}
				cat("\nRates\n")
				print(param.est)
				cat("\n")
				cat("Optima\n")
				print(theta.mat)
				cat("\n")
			}
		}
	}
	if(any(x$eigval<0)){
		index.matrix <- x$index.mat
		if(x$simmap.tree==FALSE){
			colnames(index.matrix) <- levels(x$tot.states)
		}
		if(x$simmap.tree==TRUE){
			colnames(index.matrix) <- colnames(x$phy$mapped.edge)
		}				
		#If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
		if (any(x$eigval<0)) {
			cat("The objective function may be at a saddle point -- check eigenvectors or try a simpler model", "\n")
		}
	}
	else{
		cat("Arrived at a reliable solution","\n")
	}	
}

