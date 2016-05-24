#Joint optimization of multiple traits using OUwie

#written by Jeremy M. Beaulieu

OUwie.joint <- function(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMVr","OUMA","OUMAr","OUMVA","OUMVAr"), ntraits, allfree=TRUE, simmap.tree=FALSE, scaleHeight=FALSE, root.station=TRUE, mserr="none", diagn=FALSE, quiet=FALSE){

	#Makes sure the data is in the same order as the tip labels
	if(mserr=="none" | mserr=="est"){
		data<-data.frame(data[,1], data[,2], data[,3:(2+ntraits)])
	}
	if(mserr=="known"){
		stop("You specified measurement error and it is not supported yet.")
	}
	tot.states<-factor(c(phy$node.label,as.character(data[,2])))
	k <- length(levels(tot.states))
	x <- as.matrix(data[,3:(ntraits+2)])
	if(allfree==TRUE){
		if (is.character(model)) {
			index.mat<-matrix(0,2,k*ntraits)
			if (model == "BM1"){
				np=1*ntraits
				index<-matrix(TRUE,2,k*ntraits)
				if(mserr=="est"){
					index.mat[1,1:(k*ntraits)]<-np+2
				}
				else{
					index.mat[1,1:(k*ntraits)]<-np+1				
				}
				count<-1
				for(i in seq(from = 1, by = k, length.out = ntraits)){
					j=(i+1):(i+k-1)
					index.mat[2,c(i,j)]<-count
					count<-count+1
				}
				param.count<-np+1
				bool=TRUE
			}
			#The group mean model of Thomas et al, is trivial: set bool to be TRUE:
			if (model == "BMS"){
				np=k*ntraits
				index<-matrix(TRUE,2,k*ntraits)
				if(mserr=="est"){
					index.mat[1,1:(k*ntraits)]<-np+2
				}
				else{
					index.mat[1,1:(k*ntraits)]<-np+1
				}
				index.mat[2,1:(k*ntraits)]<-1:np
				param.count<-np+1
				bool=FALSE
			}
			if (model == "OU1"){
				np=2*ntraits
				index<-matrix(TRUE,2,k*ntraits)
				for(i in seq(from = 1, by = k, length.out = ntraits)){
					j=(i+1):(i+k-1)
					index.mat[1,c(i,j)]<- max(index.mat)+1
					index.mat[2,c(i,j)]<- max(index.mat)+1
				}
				if(root.station==TRUE){
					param.count<-np+1
				}
				if(root.station==FALSE){
					param.count<-np+2
				}
				bool=root.station
			}
			if (model == "OUM"){
				np=2*ntraits
				index<-matrix(TRUE,2,k)
				for(i in seq(from = 1, by = k, length.out = ntraits)){
					j=(i+1):(i+k-1)
					index.mat[1,c(i,j)]<- max(index.mat)+1
					index.mat[2,c(i,j)]<- max(index.mat)+1
				}
				if(root.station==TRUE){
					param.count<-np+k
				}
				if(root.station==FALSE){
					param.count<-np+k+1
				}
				bool=root.station
			}
			if (model == "OUMV") {
				np=(k+1)*ntraits
				index<-matrix(TRUE,2,k)
				for(i in seq(from = 1, by = k, length.out = ntraits)){
					j=(i+1):(i+k-1)
					index.mat[1,c(i,j)] <- max(index.mat)+1
					index.mat[2,c(i,j)] <- seq(max(index.mat)+1, by=1, max(index.mat)+k)
				}	
				if(root.station==TRUE){
					param.count<-np+k
				}
				if(root.station==FALSE){
					param.count<-np+k+1
				}
				bool=root.station
			}
			if (model == "OUMVr") {
				np=(k*ntraits)+1
				index<-matrix(TRUE,2,k)
				index.mat[1,]<-1
				index.mat[2,]<-2:np				
				if(root.station==TRUE){
					param.count<-np+k
				}
				if(root.station==FALSE){
					param.count<-np+k+1
				}
				bool=root.station
			}	
			if (model == "OUMA") {
				np=(k+1)*ntraits
				index<-matrix(TRUE,2,k*ntraits)
				for(i in seq(from = 1, by = k, length.out = ntraits)){
					j=(i+1):(i+k-1)
					index.mat[1,c(i,j)]<-seq(max(index.mat)+1, by=1, max(index.mat)+k)
					index.mat[2,c(i,j)]<-max(index.mat)+1
				}
				if(root.station==TRUE){
					param.count<-np+k
				}
				if(root.station==FALSE){
					param.count<-np+k+1
				}
				bool=root.station
			}
			if (model == "OUMAr") {
				np=(k*ntraits)+1
				index<-matrix(TRUE,2,k)
				index.mat[2,]<-1
				index.mat[1,]<-2:np				
				if(root.station==TRUE){
					param.count<-np+k
				}
				if(root.station==FALSE){
					param.count<-np+k+1
				}
				bool=root.station
			}									
			if (model == "OUMVA") {
				np=k*2*ntraits
				index<-matrix(TRUE,2,k)
				index.mat[index]<-1:np
				if(root.station==TRUE){
					param.count<-np+k
				}
				if(root.station==FALSE){
					param.count<-np+k+1
				}
				bool=root.station
			}
		}
	}else{
		if (is.character(model)) {
			index.mat<-matrix(0,2,k)
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
		}
		index.mat.tmp<-index.mat
		for(i in 2:ntraits){
			index.mat<-cbind(index.mat,index.mat.tmp)
		}		
	}
	Rate.mat <- matrix(1, 2, k*ntraits)
	optim.dev <- function(p, phy, data, model, ntraits, index.mat, mserr){
		p = exp(p)
		if(model=="OUMVr" | model == "OUMAr"){
			if(model=="OUMVr"){	
				model.tmp = "OUMV"
			}
			if(model=="OUMAr"){
				model.tmp = "OUMA"
			}
		}else{
			model.tmp = model
		}
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		loglik <- 0
		count <- 3
		for(i in seq(from = 1, by = k, length.out = ntraits)){
			j=(i+1):(i+k-1)
			tmp<-NA
			try(tmp <- OUwie.fixed(phy,data[,c(1,2,count)], model=model.tmp, simmap.tree=simmap.tree, alpha=c(Rate.mat[1,c(i,j)]), sigma.sq=c(Rate.mat[2,c(i,j)]), quiet=TRUE)$loglik, silent=TRUE)
			if(!is.finite(tmp)){
				return(10000000)
			}
			if(is.na(tmp)){
				return(10000000)
			}else{
				loglik <- loglik + (-1*tmp) 
			}
			count <- count + 1
		}
		return(loglik)
	}
	
	if(quiet==FALSE){
		cat("Initializing...", "\n")
	}
	
	lb = -20
	ub = 20
	lower = rep(lb, np)
	upper = rep(ub, np)
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	
	if(model == "OU1" | model == "OUM" | model == "OUMV" | model == "OUMVr" | model == "OUMA" | model == "OUMAr" | model == "OUMVA"){
		n=max(phy$edge[,1])
		C.mat<-vcv.phylo(phy)
		a<-as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		A<-matrix(rep(a,nrow(x)),nrow(x),ncol(x),byrow=T)
		sig <- t(x-A)%*%pseudoinverse(C.mat)%*%(x-A)/(nrow(C.mat)-1)
		init.np=2*ntraits
		init.lower = rep(lb, init.np)
		init.upper = rep(ub, init.np)
		init.index.mat <- matrix(0,2,k)
		init.index.mat[1,1:k] <- 1
		init.index.mat[2,1:k] <- 2
		if(allfree==TRUE){
			init.index.mat.tmp <- init.index.mat
			#start.vals<-rep(diag(sig)[1],k)
			start.vals<-c(log(2)/max(branching.times(phy)),diag(sig)[1])
			for(i in 2:ntraits){
				init.index.mat <- cbind(init.index.mat, init.index.mat.tmp+2)
				init.index.mat.tmp <- init.index.mat.tmp+2
				#start.vals <- c(start.vals, rep(diag(sig)[i],k))
				start.vals <- c(start.vals, c(log(2)/max(branching.times(phy)),diag(sig)[i]))
			}
			init <- nloptr(x0=log(start.vals), eval_f=optim.dev, lb=init.lower, ub=init.upper, opts=opts, phy=phy, data=data, model="OU1", ntraits=ntraits, index.mat=init.index.mat, mserr="none")
			init.ip <- exp(init$solution)
			if(model == "OU1" | model == "OUM"){
				ip = init.ip
			}else{
				if(model == "OUMVr" | model == "OUMAr"){
					if(model=="OUMVr"){
						count <- 1
						ip <- mean(init.ip[seq(from = 1, by = 2, length.out = ntraits)])
						for(i in seq(from = 1, by = k, length.out = ntraits)){
							j=(i+1):(i+k-1)
							ip<-c(ip,rep(init.ip[count+1],length(unique(index.mat[2,c(i,j)]))))
							count <- count + 2
						}
					}
					if(model == "OUMAr"){
						ip<-c()
						count <- 1
						ip <- mean(init.ip[seq(from = 2, by = 2, length.out = ntraits)])
						for(i in seq(from = 1, by = k, length.out = ntraits)){
							j=(i+1):(i+k-1)
							ip <- c(ip,rep(init.ip[count],length(unique(index.mat[1,c(i,j)]))))
							count <- count + 2
						}
					}
				}else{
					ip<-c()
					if(model == "OUMA" | model == "OUMV"){ 
						count <- 1
						for(i in seq(from = 1, by = k, length.out = ntraits)){
							j=(i+1):(i+k-1)
							ip<-c(ip,rep(init.ip[count],length(unique(index.mat[1,c(i,j)]))),rep(init.ip[count+1],length(unique(index.mat[2,c(i,j)]))))
							count <- count + 2 
						}
					}
					if(model == "OUMVA"){
						count <- 1
						for(i in seq(from = 1, by = k, length.out = ntraits)){
							j=(i+1):(i+k-1)
							ip<-c(ip,rep(c(init.ip[count],init.ip[count+1]),k))
							count <- count + 2
						}
					}
				}
			}
		}
		else{
			init.np=2
			init.lower = rep(lb, init.np)
			init.upper = rep(ub, init.np)
			init.index.mat <- matrix(0,2,k)
			init.index.mat[1,1:k] <- 1
			init.index.mat[2,1:k] <- 2			
			if(model=="OU1" | model == "OUM"){
				#start.vals <- rep(mean(diag(sig)),2)
				start.vals<-c(log(2)/max(branching.times(phy)),mean(diag(sig)))
				init <- nloptr(x0=log(start.vals), eval_f=optim.dev, lb=lower, ub=upper, opts=opts, phy=phy, data=data, model="OU1", ntraits=ntraits, index.mat=init.index.mat, mserr="none")
				init.ip <- exp(c(init$solution[1],init$solution[2]))				
				ip=init.ip
			}
			else{
				#start.vals <- rep(mean(diag(sig)),np)
				start.vals <- c(log(2)/max(branching.times(phy)),mean(diag(sig)))
				init <- nloptr(x0=log(start.vals), eval_f=optim.dev, lb=init.lower, ub=init.upper, opts=opts, phy=phy, data=data, model="OU1", ntraits=ntraits, index.mat=init.index.mat, mserr="none")
				init.ip <- exp(init$solution)
				if(model == "OUMA" | model == "OUMV"){ 
					ip<-c(rep(init.ip[1],length(unique(index.mat[1,]))),rep(init.ip[2],length(unique(index.mat[2,]))))
				}
				if(model == "OUMVA"){
					ip<-rep(c(init.ip),k)
				}
			}		
		}
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		out = nloptr(x0=log(ip), eval_f=optim.dev, lb=lower, ub=upper, opts=opts, phy=phy, data=data, model=model, ntraits=ntraits, index.mat=index.mat, mserr="none")
	}
	else{
		#Starting values follow from phytools:
		C.mat<-vcv.phylo(phy)
		a<-as.numeric(colSums(solve(C.mat))%*%x/sum(solve(C.mat)))
		A<-matrix(rep(a,nrow(x)),nrow(x),ncol(x),byrow=T)
		sig <- t(x-A)%*%pseudoinverse(C.mat)%*%(x-A)/(nrow(C.mat)-1)
		#####################
		if(allfree==TRUE){
			if(model=="BM1"){
				ip=diag(sig)
			}
			if(model=="BMS"){
				ip=c()
				sigs<-diag(sig)
				for(i in 1:length(sigs)){
					ip <- c(ip, c(rep(sigs[i],k)))
				}
			}
			if(mserr=="est"){
				ip<-c(ip,0)
				lower = c(lower,0)
				upper = c(upper,10)
			}
			if(quiet==FALSE){
				cat("Finished. Begin thorough search...", "\n")
			}
			out = nloptr(x0=log(ip), eval_f=optim.dev, lb=lower, ub=upper, opts=opts, phy=phy, data=data, model=model, ntraits=ntraits, index.mat=index.mat, mserr="none")
		}else{
			if(model=="BM1"){
				ip=mean(diag(sig))
			}
			if(model=="BMS"){
				ip=c()
				for(i in 1:length(diag(sig))){
					ip <- rep(mean(diag(sig)), k)
				}
			}
			if(quiet==FALSE){
				cat("Finished. Begin thorough search...", "\n")
			}
			out = nloptr(x0=log(ip), eval_f=optim.dev, lb=lower, ub=upper, opts=opts, phy=phy, data=data, model=model, ntraits=ntraits, index.mat=index.mat, mserr="none")
		}
	}
	loglik <- -out$objective
	out$solution = exp(out$solution)
	solution<-matrix(out$solution[index.mat], dim(index.mat))
	rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq")
	if(simmap.tree==FALSE){
		colnames(solution) <- rep(levels(tot.states),ntraits)
	}
	if(simmap.tree==TRUE){
		colnames(solution) <- rep(colnames(phy$mapped.edge),ntraits)
	}
	thetas <- c()
	count <- 3
	for(i in seq(from = 1, by = k, length.out = ntraits)){
		if(model=="OUMVr" | model == "OUMAr"){
			if(model=="OUMVr"){	
				model.tmp = "OUMV"
			}
			if(model=="OUMAr"){
				model.tmp = "OUMA"
			}
		}else{
			model.tmp = model
		}
		if(model.tmp == "BM1" | model.tmp == "BMS" | model.tmp == "OU1"){
			j=(i+1):(i+k-1)
			tmp <- OUwie.fixed(phy,data[,c(1,2,count)], model=model.tmp, simmap.tree=simmap.tree, alpha=c(solution[1,c(i,j)]), sigma.sq=c(solution[2,c(i,j)]), quiet=TRUE)$theta
			tmp <- t(tmp)
			tmp.mat <- matrix(0,2,k)
			for(k in 1:k){
				tmp.mat[,k] <- tmp[,1]
			}
			thetas <- cbind(thetas,tmp.mat)
			count <- count + 1			
		}else{
			j=(i+1):(i+k-1)
			tmp <- OUwie.fixed(phy,data[,c(1,2,count)],model=model.tmp, simmap.tree=simmap.tree, alpha=c(solution[1,c(i,j)]), sigma.sq=c(solution[2,c(i,j)]), quiet=TRUE)$theta
			thetas <- cbind(thetas,t(tmp)) 
			count <- count + 1
		}
	}
	ntips = Ntip(phy)
	obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))),model=model,solution=solution, thetas=thetas, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, opts=opts, data=data, phy=phy, root.station=root.station, lb=lower, ub=upper, iterations=out$iterations, ntraits=ntraits) 
	class(obj)<-"OUwie.joint"		
	return(obj)
}


print.OUwie.joint <- function(x, ...){
	ntips = Ntip(x$phy)
	output <- data.frame(x$loglik,x$AIC,x$AICc,x$model, ntips, row.names="")
	names(output) <- c("-lnL","AIC","AICc","model","ntax")
	cat("\nOverall Fit\n")
	print(output)
	cat("\n")
	param.est <- x$solution
	rownames(x$thetas) <- c("estimate", "se")
	if(x$simmap.tree==FALSE){
		colnames(x$thetas) <- rep(levels(x$tot.states), x$ntraits)
	}
	if(x$simmap.tree==TRUE){
		colnames(x$thetas) <- rep(c(colnames(x$phy$mapped.edge)), x$ntraits)
	}
	k = length(levels(x$tot.states))
	count=1
	for(i in seq(from = 1, by = k, length.out = x$ntraits)){
		j=(i+1):(i+k-1)
		cat(paste("Trait", count, sep=" "),"\n")
		cat("\nRates\n")
		print(param.est[,c(i,j)])
		cat("\n")		
		cat("Optima\n")
		print(x$thetas[,c(i,j)])
		cat("\n")
		count<-count+1
	}
	cat("\n")
}


######################################################################################################################################
######################################################################################################################################
### Below are the model sets 
######################################################################################################################################
######################################################################################################################################

#res<-matrix(,15,20)
#res1<-MultiTraitOUwie(tree, trait, model=c("BM1"), ntraits=3, allfree=TRUE)
#res[1,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$thetas[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("BMS"), ntraits=3,allfree=TRUE)
#res[2,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$thetas[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OU1"),ntraits=3, allfree=TRUE)
#res[3,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$thetas[1,])
#write.table(res, file="fleshy.multi.set", quote=FALSE, row.names=FALSE, sep="\t")
#res1<-MultiTraitOUwie(tree, trait, model=c("OUM"), ntraits=3,allfree=TRUE)
#res[4,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMV"),ntraits=3, allfree=TRUE)
#res[5,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMVr"),ntraits=3, allfree=TRUE)
#res[6,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#write.table(res, file="fleshy.multi.set", quote=FALSE, row.names=FALSE, sep="\t")
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMA"), ntraits=3,allfree=TRUE)
#res[7,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMAr"),ntraits=3, allfree=TRUE)
#res[8,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMVA"),ntraits=3, allfree=TRUE)
#res[9,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#write.table(res, file="fleshy.multi.set", quote=FALSE, row.names=FALSE, sep="\t")
#res1<-MultiTraitOUwie(tree, trait, model=c("BMS"), ntraits=3,allfree=FALSE)
#res[10,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$thetas[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OU1"),ntraits=3, allfree=FALSE)
#res[11,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$thetas[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUM"), ntraits=3,allfree=FALSE)
#res[12,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#write.table(res, file="fleshy.multi.set", quote=FALSE, row.names=FALSE, sep="\t")
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMV"), ntraits=3,allfree=FALSE)
#res[13,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMA"), ntraits=3,allfree=FALSE)
#res[14,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#res1<-MultiTraitOUwie(tree, trait, model=c("OUMVA"), ntraits=3,allfree=FALSE)
#res[15,]<-c(res1$loglik,res1$AICc,res1$solution[1,],res1$solution[2,],res1$theta[1,])
#write.table(res, file="fleshy.multi.set", quote=FALSE, row.names=FALSE, sep="\t")



