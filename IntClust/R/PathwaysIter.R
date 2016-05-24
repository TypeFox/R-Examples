PathwaysIter<-function(List,Selection=NULL,GeneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),GeneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,niter=10,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	PathwaysOutput = list() 
	for (i in 1:niter){
		message(paste("Iteration",i,sep=" "))
		mlp = Pathways(List,Selection,GeneExpr,nrclusters,method,GeneInfo,geneSetSource,topP,topG,GENESET,sign=0.05,fusionsLog,WeightClust,names)
		PathwaysOutput [[length(PathwaysOutput )+1]] = mlp
		names(PathwaysOutput )[i]=paste("Iteration",i,sep=" ")
	 }
		
	
	return(PathwaysOutput)
}
