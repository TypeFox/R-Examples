###function to handle fixed & random (bootstrap) effects for epCA	
epCA.inference.battery <- function(DATA, DESIGN = NULL, make_design_nominal = TRUE, masses = NULL, weights = NULL, hellinger = FALSE, symmetric = TRUE, graphs = TRUE, k = 0, test.iters=100, critical.value=2){

####private functions
permute.components.ca <- function(DATA,hellinger=FALSE,symmetric=TRUE,masses=NULL,weights=NULL,k=0){
	perm.DATA <- contingency.data.break(DATA,boot=FALSE)
	return(epCA(perm.DATA,hellinger=hellinger,symmetric=symmetric,graphs=FALSE,k=k,masses=masses,weights=weights)$ExPosition.Data$eigs)
}
####end private

	fixed.res <- epCA(DATA, DESIGN, make_design_nominal, masses, weights, hellinger, symmetric, graphs=FALSE, k)

	ncomps <- fixed.res$ExPosition.Data$pdq$ng
	fj.boot.array <- array(0,dim=c(ncol(DATA), ncomps,test.iters))
	eigs.perm.matrix <- matrix(0,test.iters, min(dim(DATA)))
	
	pb <- txtProgressBar(1,test.iters,1,style=1)
	for(i in 1:test.iters){
		if(i==1){
			start.time <- proc.time()
		}
		
		fj.boot.array[,,i] <- boot.compute.fj(DATA,fixed.res,DESIGN)
		perm.eigs <- permute.components.ca(DATA,hellinger=hellinger,symmetric=symmetric,masses=masses,weights=weights,k=k)
		eigs.perm.matrix[i,1:length(perm.eigs)] <- perm.eigs
		
		if(i==1){
			cycle.time <- (proc.time() - start.time) #this is in seconds...
			if(!continueResampling(cycle.time[1] * test.iters)){
				##exit strategy.
				return(fixed.res)
			}
		}
		setTxtProgressBar(pb,i)
	}
	
#	rownames(fj.boot.array) <- colnames(DATA)
#	fj.boot.data <- list(fj.tests=boot.ratio.test(fj.boot.array,critical.value=critical.value),fj.boots=fj.boot.array)
#	class(fj.boot.data) <- c("inpoBootstrap", "list")
	rownames(fj.boot.array) <- colnames(DATA)
	boot.ratio.test.data <- boot.ratio.test(fj.boot.array,critical.value=critical.value)
	class(boot.ratio.test.data) <- c("inpoBootTests","list")
	fj.boot.data <- list(tests=boot.ratio.test.data,boots=fj.boot.array)
	class(fj.boot.data) <- c("inpoBoot", "list")
	
	##do I still need this rounding?
	eigs.perm.matrix <- round(eigs.perm.matrix,digits=15)	
	inertia.perm <- rowSums(eigs.perm.matrix)
	omni.p <- max(1-(sum(inertia.perm < sum(fixed.res$ExPosition.Data$eigs))/test.iters),1/test.iters)
	omni.data <- list(p.val=round(omni.p,digits=4),inertia.perm=inertia.perm, inertia=sum(fixed.res$ExPosition.Data$eigs))
	class(omni.data) <- c("inpoOmni","list")
		
	eigs.perm.matrix <- eigs.perm.matrix[,1:ncomps]
	component.p.vals <- 1-(colSums(eigs.perm.matrix < matrix(fixed.res$ExPosition.Data$eigs,test.iters, ncomps,byrow=TRUE))/test.iters)
	component.p.vals[which(component.p.vals==0)] <- 1/test.iters	
	components.data <- list(p.vals=round(component.p.vals,digits=4), eigs.perm=eigs.perm.matrix, eigs=fixed.res$ExPosition.Data$eigs)
	class(components.data) <- c("inpoComponents","list")
				
 	Inference.Data <- list(components=components.data,fj.boots=fj.boot.data,omni.data=omni.data)
	class(Inference.Data) <- c("epCA.inference.battery","list")

	ret.data <- list(Fixed.Data=fixed.res,Inference.Data=Inference.Data)
	class(ret.data) <- c("inpoOutput","list")	

	#graphing needs to happen here.
	if(graphs){
		inGraphs(ret.data)
	}

	return(ret.data)
	
}