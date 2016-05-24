PathwayAnalysis<-function(List,Selection=NULL,GeneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),GeneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,niter=10,fusionsLog=TRUE,WeightClust=TRUE,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
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
	
	Pathways=PathwaysIter(List,Selection,GeneExpr,nrclusters,method,GeneInfo,geneSetSource,topP,topG,GENESET,sign,niter,fusionsLog,WeightClust,names)
	
	if(is.null(Selection)){
		Selection=FALSE
	}
	else{
		Selection=TRUE
	}
	
	Intersection=Geneset.intersect(PathwaysOutput=Pathways,Selection,sign,names,seperatetables,separatepvals)
	
	return(Intersection)
	
}