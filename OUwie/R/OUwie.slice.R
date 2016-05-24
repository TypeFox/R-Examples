#OUwie when selective regimes are time slices

#written by Jeremy M. Beaulieu

OUwie.slice<-function(phy, data, model=c("BMS","OUM","OUMV","OUMA","OUMVA"), timeslices=c(NA), scaleHeight=FALSE, root.station=TRUE, mserr="none", slice.lower.bound=NULL, diagn=FALSE, quiet=FALSE, warn=TRUE){
	
	#Makes sure the data is in the same order as the tip labels
	if(mserr=="none" | mserr=="est"){
		data<-data.frame(data[,2], data[,2], row.names=data[,1])
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
	simmap.tree<-TRUE
	#How many regimes and timeslices from the input
	if(any(timeslices>0, na.rm=T) | any(is.na(timeslices))){
		k<-length(timeslices)+1
	}
	else{
		stop("You have not specified a timeslice")
	}
	max.height <- max(nodeHeights(phy))
	timeslices <- max.height - timeslices
	timeslices <- c(0,timeslices)
	#Values to be used throughout
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	
	#A boolean for whether the root theta should be estimated -- default is that it should be.
	root.station=root.station
	#Obtains the state at the root
	root.state=1
	##Begins the construction of the edges matrix -- similar to the ouch format##
	edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
	
	if(scaleHeight==TRUE){
		edges[,4:5]<-edges[,4:5]/max.height
	}
	edges=edges[sort.list(edges[,3]),]
	
	#Resort the edge matrix so that it looks like the original matrix order
	edges=edges[sort.list(edges[,1]),]
	x<-as.matrix(data[,2])
	#Matches the model with the appropriate parameter matrix structure
	if (is.character(model)) {
		index.mat<-matrix(0,2,k)
		if (model == "BMS"){
			np=k
			index<-matrix(TRUE,2,k)
			if(mserr=="est"){
				index.mat[1,1:k]<-np+2
			}
			else{
				index.mat[1,1:k]<-0
			}
			index.mat[2,1:k]<-1:np
			if(root.station==TRUE){
				param.count<-np+k
			}
			if(root.station==FALSE){
				param.count<-np+1
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
	}
	if(warn==TRUE){
		if(param.count > (ntips/10)){
			warning("You might not have enough data to fit this model well", call.=FALSE, immediate.=TRUE)
		}
	}
	Rate.mat <- matrix(1, 2, k)
	#Generates vector for identifying any timeslices that need to be estimated
	Slices.vector <- timeslices
	if (length(which(is.na(timeslices)))>0) {
		Slices.vector[which(is.na(timeslices))] <- max(index.mat)+sequence(length(which(is.na(timeslices))))
		Slices.vector[which(!is.na(timeslices))] <- 0
		index.mat[index.mat==0] <- Slices.vector[Slices.vector == 0] <- max(Slices.vector)+sequence(length(which(is.na(timeslices))))
	}
	else{
		Slices.vector[which(!is.na(timeslices))] <- 0
		index.mat[index.mat==0] <- Slices.vector[Slices.vector == 0] <- max(index.mat)+1
	}
	#Likelihood function for estimating model parameters
	dev.slice<-function(p, index.mat, timeslices, mserr){
		p = exp(p)
		non.estimated.slices<-timeslices[which(!is.na(timeslices))]
		timeslices[] <-c(p, 0)[Slices.vector]
		timeslices[timeslices==0]<-non.estimated.slices
		timeslices<-timeslices[order(timeslices, decreasing=FALSE)]
		if(scaleHeight==TRUE){
			timeslices = timeslices/max.height
		}
		phy.sliced<-make.era.map(phy,timeslices)
		Rate.mat[] <- c(p, 1e-12)[index.mat]
		N<-length(x[,1])
		V<-varcov.ou(phy.sliced, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight)
		W<-weight.mat(phy.sliced, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, assume.station=bool)
		if (any(is.nan(diag(V))) || any(is.infinite(diag(V)))) return(10000000)		
		if(mserr=="known"){
			diag(V)<-diag(V)+data[,3]
		}
		if(mserr=="est"){
			diag(V)<-diag(V)+p[length(p)]
		}
		theta<-Inf
		try(theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x, silent=TRUE)
		if(any(theta==Inf)){
			return(10000000)
		}
		#Seemingly more stable calculation of the determinant:
		DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
		#However, sometimes there are underflow issues so I toggle between the two approaches:
		if(!is.finite(DET)){
			DET<-determinant(V, logarithm=TRUE)
			logl<--.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
		}else{
			logl<--.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET)-.5*(N*log(2*pi))
		}
		if(!is.finite(logl)){
			return(10000000)
		}		
		return(-logl)
	}
	
	if(quiet==FALSE){
		cat("Initializing...", "\n")
	}
    if(is.null(slice.lower.bound)){
        slice.lower.bound = 0.00001
    }    
	lb = -20
	ub = 20
	lower = rep(lb, np)
	upper = rep(ub, np)
	#Update the bounds with estimated timeslices:
	lower = c(lower, rep(-21,length(which(is.na(timeslices)))))
	upper = c(upper, log(rep(max(nodeHeights(phy))-slice.lower.bound,length(which(is.na(timeslices))))))
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	if(model == "OUM" | model == "OUMV" | model == "OUMA" | model == "OUMVA"){
		#Need to construct a dataset so that we can run OUwie to get starting values:
		data.tmp <- data.frame(Genus_species=phy$tip.label,reg=rep(1,length(data[,1])),contT=data[,2])
		phy.tmp <- phy
		phy.tmp$node.label <- sample(c(1:k), phy.tmp$Nnode, replace=T)
		init <- OUwie(phy.tmp, data.tmp, model="OU1", simmap.tree=FALSE, scaleHeight=scaleHeight, root.station=root.station, mserr=mserr, diagn=FALSE, quiet=TRUE)
		init.ip <- c(init$solution[1,1],init$solution[2,1])

		if(model=="OUMV" | model=="OUMA" | model=="OUM"){
			ip<-c(rep(init.ip[1],length(unique(index.mat[1,]))),rep(init.ip[2],length(unique(index.mat[2,]))))
		}
		if(model=="OUMVA"){
			ip<-c(rep(init.ip,k))
		}

		if(mserr=="est"){
			ip <- c(ip,0)
			lower = c(lower,0)
			upper = c(upper,ub)
		}
		if(any(is.na(timeslices))){
			ip <- c(ip, max(nodeHeights(phy))/2)
		}
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		out = nloptr(x0=log(ip), eval_f=dev.slice, lb=lower, ub=upper, opts=opts, index.mat=index.mat, timeslices=timeslices, mserr=mserr)
	}
	else{
		if(scaleHeight==TRUE){
			d <- max(diag(vcv.phylo(phy)))
			phy$edge.length<-(phy$edge.length/d)
		}		
		#Starting values follow from phytools:
		C.mat<-vcv.phylo(phy)
		a<-as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		A<-matrix(rep(a,nrow(x)),nrow(x),ncol(x), byrow=TRUE)
		sig<-as.numeric(t(x-A)%*%solve(C.mat)%*%(x-A)/n)
		ip=rep(sig,k)
		if(mserr=="est"){
			ip<-c(ip,0)
			lower = c(lower,0)
			upper = c(upper,10)
		}
		if(any(is.na(timeslices))){
			ip <- c(ip, max(nodeHeights(phy))/2)
		}
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		out = nloptr(x0=log(ip), eval_f=dev.slice, lb=lower, ub=upper, opts=opts, index.mat=index.mat, timeslices=timeslices, mserr=mserr)
	}
	
	loglik <- -out$objective
	out$solution = exp(out$solution)
	#Takes estimated parameters from dev and calculates theta for each regime:
	dev.theta.slice<-function(p, index.mat, timeslices=timeslices, mserr=mserr){
		non.estimated.slices<-timeslices[which(!is.na(timeslices))]
		timeslices[] <- c(p, 0)[Slices.vector]
		timeslices[timeslices==0] <- non.estimated.slices
		timeslices <- timeslices[order(timeslices, decreasing=FALSE)]
		phy.sliced <- make.era.map(phy,timeslices)
		tmp<-NULL
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		N<-length(x[,1])
		V<-varcov.ou(phy.sliced, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight)
		W<-weight.mat(phy.sliced, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, assume.station=bool)
		
		if(mserr=="known"){
			diag(V)<-diag(V)+data[,3]
		}
		if(mserr=="est"){
			diag(V)<-diag(V)+p[length(p)]
		}
		theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
		#Calculates the hat matrix:
		#H.mat<-W%*%pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)
		#Standard error of theta -- uses pseudoinverse to overcome singularity issues
		se<-sqrt(diag(pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)))
		#tmp$res<-W%*%theta-x
		#Joins the vector of thetas with the vector of standard errors into a 2 column matrix for easy extraction at the summary stage
		tmp$theta.est<-cbind(theta,se)
		#Returns final GLS solution
		tmp
	}
	#Informs the user that the summarization has begun, output model-dependent and dependent on whether the root theta is to be estimated
	if(quiet==FALSE){
		cat("Finished. Summarizing results.", "\n")	
	}
	theta <- dev.theta.slice(out$solution, index.mat, timeslices, mserr)
	non.estimated.slices<-timeslices[which(!is.na(timeslices))]
	timeslices[] <- c(out$solution, 0)[Slices.vector]
	timeslices[timeslices==0] <- non.estimated.slices
	timeslices <- timeslices[order(timeslices, decreasing=FALSE)]
	phy.sliced<-make.era.map(phy,timeslices)
	tot.states<-factor(colnames(phy.sliced$mapped.edge))
	#Calculates the Hessian for use in calculating standard errors and whether the maximum likelihood solution was found
	if(diagn==TRUE){
		h <- hessian(x=out$solution, func=dev.slice, index.mat=index.mat, timeslices=timeslices, mserr=mserr)
		#Using the corpcor package here to overcome possible NAs with calculating the SE
		solution <- matrix(out$solution[index.mat], dim(index.mat))
		solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[index.mat], dim(index.mat))
		rownames(solution) <- rownames(solution.se) <- rownames(index.mat) <- c("alpha","sigma.sq")
		if(simmap.tree==FALSE){
			colnames(solution) <- colnames(solution.se) <- levels(tot.states)
		}
		if(simmap.tree==TRUE){
			colnames(solution) <- colnames(solution.se) <- c(colnames(phy.sliced$mapped.edge))
		}				
		#Eigendecomposition of the Hessian to assess reliability of likelihood estimates
		hess.eig <- eigen(h,symmetric=TRUE)
		#If eigenvect is TRUE then the eigenvector and index matrix will appear in the list of objects 
		eigval <- signif(hess.eig$values,2)
		eigvect <- round(hess.eig$vectors, 2)
		if(mserr=="est"){
			mserr.est<-out$solution[length(out$solution)]
			param.count<-param.count+1
		}
		else{
			mserr.est<-NULL
		}	
		obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))),model=model,solution=solution, theta=theta$theta.est, solution.se=solution.se, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, opts=opts, data=data, phy=phy.sliced, root.station=root.station, lb=lower, ub=upper, iterations=out$iterations, res=theta$res, eigval=eigval, eigvect=eigvect) 
	}
	if(diagn==FALSE){
		solution <- matrix(out$solution[index.mat], dim(index.mat))
		rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq")
		if(simmap.tree==FALSE){
			colnames(solution) <- levels(tot.states)
		}
		if(simmap.tree==TRUE){
			colnames(solution) <- c(colnames(phy.sliced$mapped.edge))
		}
		if(mserr=="est"){
			mserr.est<-out$solution[length(out$solution)]
			param.count<-param.count+1
		}
		else{
			mserr.est<-NULL
		}
		obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))),model=model,solution=solution, theta=theta$theta.est, timeslices=timeslices, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, opts=opts, data=data, phy=phy.sliced, root.station=root.station, lb=lower, ub=upper, iterations=out$iterations, res=theta$res) 
	}
	class(obj)<-"OUwie.slice"		
	return(obj)
}

print.OUwie.slice<-function(x, ...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$model,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","model","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	if (is.character(x$model)) {
		if (x$model == "BMS"){
			param.est <- x$solution
			if(x$root.station==FALSE){
				theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
			}
			else{
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
			}
			rownames(theta.mat)<-c("estimate", "se")
			if(x$simmap.tree==FALSE){
				colnames(theta.mat) <- levels(x$tot.states)
			}
			if(x$simmap.tree==TRUE){
				colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
			}
			cat("Time slices\n")
			time.mat<-matrix(max(nodeHeights(x$phy))-x$timeslices,1,length(x$timeslices))
			colnames(time.mat) <- c(colnames(x$phy$mapped.edge))
			rownames(time.mat) <- "interval start"
			print(time.mat)
			cat("Rates\n")
			print(param.est)
			cat("\n")
			cat("Optima\n")
			print(theta.mat)
			cat("\n")
		}
		if (x$root.station == TRUE){
			if (x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
				param.est<- x$solution
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat) <- levels(x$tot.states)
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
				}
				cat("Time slices\n")
				time.mat<-matrix(max(nodeHeights(x$phy))-x$timeslices,1,length(x$timeslices))
				colnames(time.mat) <- c(colnames(x$phy$mapped.edge))
				rownames(time.mat) <- "interval start"
				print(time.mat)
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
					colnames(theta.mat)<-c("Root", levels(x$tot.states))
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat)<-c("Root", colnames(x$phy$mapped.edge))
				}
				cat("Time slices\n")
				time.mat<-matrix(max(nodeHeights(x$phy))-x$timeslices,1,2)
				colnames(time.mat) <- c(colnames(x$phy$mapped.edge))
				rownames(time.mat) <- "interval start"
				print(time.mat)
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

