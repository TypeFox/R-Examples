epPCA.inference.battery <- function(DATA, scale = TRUE, center = TRUE, DESIGN = NULL, make_design_nominal = TRUE, graphs = TRUE, k = 0, test.iters=100, constrained=FALSE,critical.value=2){

###some private functions
permute.components.pca <- function(DATA,scale=TRUE,center=TRUE,k=0){
	perm.DATA <- apply(DATA,2,sample)
	return(epPCA(perm.DATA,graphs=FALSE,scale=scale,center=center,k=k)$ExPosition.Data$eigs)
}
###end private functions
    	
	fixed.res <- epPCA(DATA, scale, center, DESIGN, make_design_nominal, graphs=FALSE, k)

	ncomps <- fixed.res$ExPosition.Data$pdq$ng
	fj.boot.array <- array(0,dim=c(ncol(DATA),ncomps,test.iters))
	eigs.perm.matrix <- matrix(0,test.iters,min(dim(DATA)))
		
	pb <- txtProgressBar(1,test.iters,1,style=1)
	for(i in 1:test.iters){
		if(i==1){
			start.time <- proc.time()
		}
				
		fj.boot.array[,,i] <- boot.compute.fj(DATA,res=fixed.res,DESIGN=DESIGN,constrained=constrained)
		perm.eigs <- permute.components.pca(DATA,scale,center,k=k)
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

	#fj.boot.array <- replace(fj.boot.array,is.nan(fj.boot.array),0)
	rownames(fj.boot.array) <- colnames(DATA)
	boot.ratio.test.data <- boot.ratio.test(fj.boot.array,critical.value=critical.value)
	class(boot.ratio.test.data) <- c("inpoBootTests","list")
	fj.boot.data <- list(tests=boot.ratio.test.data,boots=fj.boot.array)
	class(fj.boot.data) <- c("inpoBoot", "list")
	
	eigs.perm.matrix <- eigs.perm.matrix[,1:ncomps]
	component.p.vals <- 1-(colSums(eigs.perm.matrix < matrix(fixed.res$ExPosition.Data$eigs,test.iters, ncomps,byrow=TRUE))/test.iters)
	component.p.vals[which(component.p.vals==0)] <- 1/test.iters
	components.data <- list(p.vals=round(component.p.vals,digits=4), eigs.perm=eigs.perm.matrix, eigs=fixed.res$ExPosition.Data$eigs)
	class(components.data) <- c("inpoComponents","list")
	
	Inference.Data <- list(components=components.data,fj.boots=fj.boot.data)
	class(Inference.Data) <- c("epPCA.inference.battery","list")

	ret.data <- list(Fixed.Data=fixed.res,Inference.Data=Inference.Data)
	class(ret.data) <- c("inpoOutput","list")	

	if(graphs){
		inGraphs(ret.data)
	}

	return(ret.data)
}