#JOINT OPTIMIZATION OF REGIMES AND OUwie PARAMETERS.

#written by Jeremy M. Beaulieu
#source("pathLike.R")
#library(corHMM)
#library(OUwie)
#This function is intended as a framework for simultaneously optimizing a discrete character and the continuous character when running OUwie. The idea behind this is that we assume a correlation between the two trait types, but we are not modeling them in that way. In other words, we reconstruct the regimes first, then fit OUwie separately. However, where rates change on a branch may not be optimal without information from how the regimes evolved; where regimes change may not be optimal without information from how rates change. This method should not be confused with the "hypothesis-free" approaches like SURFACE or AUTEUR. Such models are useful, though in my view, the stronger alpha is, the less likely the model will be in reconstructing regimes correctly. Moreover, with methods like Brownie and OUwie, we often go into these analyses with an a priori hypothesis about what might be governing the rates and optima. So hOUwie is intended to meet everyone halfway in this regard, as well as to fix an issue that seems to have gone unnoticed as to how these analyses are conducted (though see Revell 2013).
#The basic idea is as follows: with the joint reconstruction we can find the likeliest state for each node in the tree using the algorithm of Pupko et al (2000). However, to obtain the true likelihood of a reconstruction, we need to also assess changes along a branch as well as a node. To do so is simple: simply add in a set of nodes of degree 2 down a branch. Unfortunately, there is this strange property where change occurs either at the instant of speciation or right before -- nowhere else. Need to think of a solution for this, as I find it unclear as to whether OUwie can truly provide a pushback when optimizing corHMM. In other words, for a set of parameters the reconstructions could be very similar, which would mean the OUwie estimates are likely to be the same. The obvious fix is to include a parameter for every branch that modulates where changes occur on along branch. Obviously, this is insane. But there has to be a neat trick to deal with this.

##To do:
# 1. Do optimization in log scale 
# 2. Get rid of this diagnostic stuff.
# 3. Alternative: run a set of stochastic maps, average the likelihood, and do it that way.
# 4. 

hOUwie <- function(phy, data, ouwie.model=c("BMS","OUM","OUMV","OUMA","OUMVA"), discrete.model=c("ER", "SYM", "ARD", "HRM"), scaleHeight=FALSE, root.station=TRUE, root.p=NULL, rate.cat=NULL, ntraits=1, nstates=2, rate.mat=NULL, lower.bounds=c(1e-6,0), upper.bounds=c(1000,1000), mserr="none", diagn=FALSE, quiet=FALSE, warn=TRUE){
	
	ntips <- Ntip(phy)
	if(warn==TRUE){
		if(param.count > (ntips/10)){
			warning("You might not have enough data to fit this model well", call.=FALSE, immediate.=TRUE)
		}
	}

	if(discrete.model=="HRM"){
		hrm=TRUE
	}else{
		hrm=FALSE
	}
	
	if(is.null(rate.mat)){
		rate.mat <- rate.mat.maker(rate.cat=rate.cat, hrm=hrm, ntraits=ntraits, nstates=nstates, model=discrete.model)
		np.discrete <- max(rate.mat, na.rm=TRUE)
		lower.discrete = rep(lower.bounds[2], np.discrete)
		upper.discrete = rep(upper.bounds[2], np.discrete)
	}else{
		np.discrete <- max(rate.mat, na.rm=TRUE)
		lower.discrete = rep(lower.bounds[2], np.discrete)
		upper.discrete = rep(upper.bounds[2], np.discrete)
	}
	
	dev.joint <- function(p, phy, data, mserr, discrete.model, ouwie.model, np.discrete, np.continuous, hrm, rate.cat, ntraits, rate.mat, root.p, root.station, optimize){
		param.count <- np.discrete+np.continuous
		if(ouwie.model=="BMS"){
			alpha<-rep(1e-10,dim(rate.mat)[1])
			sigma.sq <- p[1:np.continuous]
		}
		if(ouwie.model=="OUM"){
			alpha <- rep(p[1], dim(rate.mat)[1])
			sigma.sq <- rep(p[2], dim(rate.mat)[1])
		}
		if(ouwie.model=="OUMV"){
			alpha <- rep(p[1], dim(rate.mat)[1])
			sigma.sq <- p[2:np.continuous]
		}
		if(ouwie.model=="OUMA"){
			alpha <- p[1:(np.continuous-1)]
			sigma.sq <- p[np.continuous]
		}
		if(ouwie.model=="OUMVA"){
			alpha <- p[1:(np.continuous/2)]
			sigma.sq <- p[((np.continuous/2)+1):np.continuous]
		}
		trans.rate <- p[(np.continuous+1):param.count]
		print(trans.rate)
		regime.lik <- PathLik(phy, data[,c(1:2)], p=trans.rate, hrm=hrm, rate.cat=rate.cat, ntraits=ntraits, rate.mat=rate.mat, model=discrete.model, root.p=root.p)
		print(paste("corHMM", regime.lik$loglik))
		phy.painted <- regime.lik$mapped.tree
		print(colnames(phy.painted$mapped.edge))
		print(sigma.sq[as.numeric(colnames(phy.painted$mapped.edge))])
		ouwie.lik <- OUwie.fixed(phy.painted, data, model=ouwie.model, simmap.tree=TRUE, alpha=alpha[as.numeric(colnames(phy.painted$mapped.edge))], sigma.sq=sigma.sq[as.numeric(colnames(phy.painted$mapped.edge))], root.station=root.station, mserr=mserr, quiet=TRUE)
		print(paste("ouwie", ouwie.lik$loglik))
		if(optimize==TRUE){
			loglik <- sum(regime.lik$loglik, ouwie.lik$loglik)
			print(loglik)
			return(-loglik)
		}else{
			obj<-NULL
			obj$regime.paint <- regime.lik$mapped.tree
			obj$theta <- ouwie.lik$theta
			obj$ouwie.mat <- ouwie.lik$solution
			return(obj)
		}
	}
	
	if(quiet==FALSE){
		cat("Initializing...", "\n")
	}
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	
	init.regime <- rayDISC(phy, data[,c(1,2)], model="ER", node.states="joint", root.p=root.p)

	if(mserr=="none" | mserr=="est"){
		data.sort<-data.frame(data[,2], data[,3], row.names=data[,1])
		data.sort<-data.sort[phy$tip.label,]
	}
	if(mserr=="known"){
		if(!dim(data)[2]==4){
			stop("You specified measurement error should be incorporated, but this information is missing")
		}
		else{
			data.sort<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
			data.sort<-data.sort[phy$tip.label,]
		}
	}
	x<-data[,3]
	n=max(phy$edge[,1])	
	if(ouwie.model == "OUM" | ouwie.model == "OUMV" | ouwie.model == "OUMA" | ouwie.model == "OUMVA"){
		if(scaleHeight == TRUE){
			d <- max(diag(vcv.phylo(phy)))
			phy$edge.length<-(phy$edge.length/d)
		}
		C.mat<-vcv.phylo(phy)
		a<-as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		sig<-as.numeric(t(x-a)%*%solve(C.mat)%*%(x-a)/n)
		num.states <- dim(rate.mat)[1]
		phy.tmp <- phy
		phy.tmp$node.label <- rep(c(1:num.states), Nnode(phy.tmp))
		trait.tmp <- data
		trait.tmp[,2] <- rep(1, Ntip(phy.tmp))
		init <- OUwie(phy.tmp, trait.tmp, model="OU1", quiet=TRUE, warn=FALSE)
		init.ip <- c(init$solution[1], init$solution[2])
		if(ouwie.model=="OUM"){
			ip = init.ip
			np.continuous=length(ip)
			lower.continuous = rep(lower.bounds[1], np.continuous)
			upper.continuous = rep(upper.bounds[1], np.continuous)
		}
		if(ouwie.model=="OUMV"){
			ip<-c(init.ip[1],rep(init.ip[2],dim(rate.mat)[1]))
			np.continuous=length(ip)
			lower.continuous = rep(lower.bounds[1], np.continuous)
			upper.continuous = rep(upper.bounds[1], np.continuous)
		}
		if(ouwie.model=="OUMA"){
			ip<-c(rep(init.ip[1],dim(rate.mat)[1]),init.ip[2])
			np.continuous=length(ip)
			lower.continuous = rep(lower.bounds[1], np.continuous)
			upper.continuous = rep(upper.bounds[1], np.continuous)
		}
		if(ouwie.model=="OUMVA"){
			ip<-c(rep(init.ip[1],dim(rate.mat)[1]),rep(init.ip[2],dim(rate.mat)[1]))
			np.continuous=length(ip)
			lower.continuous = rep(lower.bounds[1], np.continuous)
			upper.continuous = rep(upper.bounds[1], np.continuous)
		}
		
		lower.final <- c(lower.continuous, lower.discrete)
		upper.final <- c(upper.continuous, upper.discrete)
		
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		ip.final <- c(ip, rep(init.regime$solution[2], np.discrete))
		
		#NOTE: This might be too much stuff to pass (I cannot remember if there is an upper limit or not) -- if failure, then create object and put this stuff in it like deep42:
		out = nloptr(x0=ip.final, eval_f=dev.joint, lb=lower.final, ub=upper.final, opts=opts, phy=phy, data=data, mserr=mserr, discrete.model=discrete.model, ouwie.model=ouwie.model, np.discrete=np.discrete, np.continuous=np.continuous, hrm=hrm, rate.cat=rate.cat, ntraits=ntraits, rate.mat=rate.mat, root.p=root.p, root.station=root.station, optimize=TRUE)
	}else{
		if(scaleHeight==TRUE){
			d <- max(diag(vcv.phylo(phy)))
			phy$edge.length<-(phy$edge.length/d)
		}
		##Starting values follow directly from phytools:
		C.mat <- vcv.phylo(phy)
		a <- as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		sig <- as.numeric(t(x-a)%*%solve(C.mat)%*%(x-a)/n)
		#####################
		ip=rep(sig,dim(rate.mat)[1])
		np.continuous=length(ip)
		lower.continuous = rep(lower.bounds[1], np.continuous)
		upper.continuous = rep(upper.bounds[1], np.continuous)			
		
		lower.final <- c(lower.continuous, lower.discrete)
		upper.final <- c(upper.continuous, upper.discrete)		
		
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		ip.final <- c(ip, rep(init.regime$solution[2], np.discrete))
		out = nloptr(x0=ip.final, eval_f=dev.joint, lb=lower.final, ub=upper.final, opts=opts, phy=phy, data=data, mserr=mserr, discrete.model=discrete.model, ouwie.model=ouwie.model, np.discrete=np.discrete, np.continuous=np.continuous, hrm=hrm, rate.cat=rate.cat, ntraits=ntraits, rate.mat=rate.mat, root.p=root.p, root.station=root.station, optimize=TRUE)
	}
	
	#Optimize states and estimate proper thetas using MLEs:
	maps.and.optima <- dev.joint(p=out$solution, phy=phy, data=data, mserr=mserr, discrete.model=discrete.model, ouwie.model=ouwie.model, np.discrete=np.discrete, np.continuous=np.continuous, hrm=hrm, rate.cat=rate.cat, ntraits=ntraits, rate.mat=rate.mat, root.p=root.p, root.station=root.station, optimize=FALSE)
	param.count <- np.continuous+np.discrete

	#Formats the discrete rate matrix properly -- the OUwie rate matrix should be returned from OUwie already formatted correctly (go me!):
	trans.par <- out$solution[(np.continuous+1):param.count]
	rate.mat[is.na(rate.mat)] <- max(rate.mat,na.rm=TRUE)+1
	discrete.solution <- matrix(trans.par[rate.mat], dim(rate.mat))
	rownames(discrete.solution) <- colnames(discrete.solution) <- rownames(rate.mat)
	loglik = -out$objective
	print(maps.and.optima$ouwie.mat)
	#Generates the hOUwie object -- be sure that you include as much as needs to be there:
	obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(Ntip(phy)/(Ntip(phy)-param.count-1))), ouwie.model=ouwie.model, discrete.model=discrete.model, solution=out$solution, mserr=mserr, theta=maps.and.optima$theta, tot.states=dim(rate.mat)[1], discrete.mat=discrete.solution, ouwie.mat=maps.and.optima$ouwie.mat, opts=opts, data=data, phy=maps.and.optima$regime.paint, root.station=root.station, lb=lower.final, ub=upper.final, iterations=out$iterations) 
	class(obj)<- "hOUwie"
	return(obj)
}


print.hOUwie<-function(x, ...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik, x$AIC, x$AICc, x$discrete.model, x$ouwie.model, row.names="")
	names(output)<-c("-lnL","AIC","AICc","disc.model","ouwie.model")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	cat("Discrete Rates\n")
	print(x$discrete.mat)
	cat("\n")

	if (is.character(x$ouwie.model)) {
		if (x$ouwie.model == "BMS"){
			param.est <- x$ouwie.mat
			rownames(param.est)<-c("alpha","sigma.sq")
			theta.mat <- matrix(t(x$theta[1,]), 2, x$tot.states)
			rownames(theta.mat)<-c("estimate", "se")
			regime.names <- as.character(colnames(x$discrete.mat))
			colnames(param.est) <- colnames(theta.mat) <- regime.names[as.numeric(colnames(x$phy$mapped.edge))]
			cat("OUwie Rates\n")
			print(param.est)
			cat("\n")
			cat("Optima\n")
			print(theta.mat)
			cat("\n")
		}
		if (x$root.station == TRUE){
			if (x$ouwie.model == "OUM"| x$ouwie.model == "OUMV"| x$ouwie.model == "OUMA" | x$ouwie.model == "OUMVA"){
				param.est<- x$ouwie.mat
				rownames(param.est)<-c("alpha","sigma.sq")
				theta.mat<-matrix(t(x$theta), 2, x$tot.states)
				rownames(theta.mat)<-c("estimate", "se")
				regime.names <- colnames(x$discrete.mat)
				colnames(param.est) <- colnames(theta.mat) <- regime.names[as.numeric(colnames(x$phy$mapped.edge))]
				cat("OUwie Rates\n")
				print(param.est)
				cat("\n")
				cat("Optima\n")
				print(theta.mat)
				cat("\n")
			}
		}
		if (x$root.station == FALSE){
			if (x$ouwie.model == "OUM"| x$ouwie.model == "OUMV"| x$ouwie.model == "OUMA" | x$ouwie.model == "OUMVA"){ 
				param.est<- x$ouwie.mat
				rownames(param.est)<-c("alpha","sigma.sq")
				theta.mat<-matrix(t(x$theta), 2, x$tot.states+1)
				rownames(theta.mat)<-c("estimate", "se")
				regime.names <- colnames(x$discrete.mat)
				colnames(param.est) <- colnames(theta.mat) <- regime.names[as.numeric(colnames(x$phy$mapped.edge))]
				colnames(theta.mat)<-c("Root", regime.names[as.numeric(colnames(x$phy$mapped.edge))])
				cat("OUwie Rates\n")
				print(param.est)
				cat("\n")
				cat("Optima\n")
				print(theta.mat)
				cat("\n")
			}
		}		
	}		
}


#TO DO: Make a quiet option for corHMM. 
#Need to change namespace stuff in corHMM to allow for rayDISC functions to be seen by hOUwie.
#Also need to generate a dynamic way of figuring out which data column is continuous trait, and which one(s) are the discrete.
#Just had a thought: if the root switches states, it will not be consistent with the current par estimates.
#Fix rate.mat.maker so that it outputs colnames and rownames when characters take on states.

