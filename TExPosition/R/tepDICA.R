#tepDICA <- function(DATA,make_data_nominal=FALSE,DESIGN=NULL,make_design_nominal=TRUE,group.masses=NULL,ind.masses=NULL,weights=NULL,hellinger=FALSE,symmetric=TRUE,graphs=TRUE,k=0){
#tepDICA <- function(DATA,make_data_nominal=FALSE,DESIGN=NULL,make_design_nominal=TRUE,group.masses=NULL,ind.masses=NULL,weights=NULL,symmetric=TRUE,graphs=TRUE,k=0){
	
tepDICA <- function(DATA,make_data_nominal=FALSE,DESIGN=NULL,make_design_nominal=TRUE,group.masses=NULL,weights=NULL,symmetric=TRUE,graphs=TRUE,k=0){
	
		
	DESIGN <- texpoDesignCheck(DATA,DESIGN,make_design_nominal,force_bary=TRUE)	
	main <- deparse(substitute(DATA))	
	DATA <- as.matrix(DATA)
	if(make_data_nominal){
		DATA <- makeNominalData(DATA)
	}
	
	#Make the group x variable contingency table
	R <- t(DESIGN) %*% DATA
	colnames(R) <- colnames(DATA)
	rownames(R) <- colnames(DESIGN)	
	Rdesign <- diag(nrow(R))
	rownames(Rdesign) <- rownames(R)

	#The results from the group x variable matrix
	#res <- coreCA(R,masses=group.masses,weights=weights,hellinger=hellinger,symmetric=symmetric,k=k)
	#res <- epCA(R, DESIGN=Rdesign, make_design_nominal=FALSE, masses = group.masses, weights = weights, hellinger = hellinger, symmetric = symmetric, graphs = FALSE,k=k)
	res <- epCA(R, DESIGN=Rdesign, make_design_nominal=FALSE, masses = group.masses, weights = weights, hellinger = FALSE, symmetric = symmetric, graphs = FALSE,k=k)	
	res <- res$ExPosition.Data

	supplementaryRes <- supplementaryRows(DATA,res)
	res$fii <- supplementaryRes$fii
	res$dii <- supplementaryRes$dii
	res$rii <- supplementaryRes$rii		
	
	res$ly <- supplementaryCols(t(DESIGN),res)$fjj * (1/sqrt(sum(R)))
	res$lx <- res$fii * matrix(rowSums(DATA),nrow(res$fii),ncol(res$fii)) * (1/sqrt(sum(R)))
		
	assignments <- fii2fi(DESIGN,res$fii,res$fi)
	assignments$r2 <- R2(res$M,res$di,ind.masses=NULL,res$dii)
	class(assignments) <- c("tepAssign","list")
	res$assign <- assignments

	#new res here
	class(res) <- c("tepDICA","list")	
	tepPlotInfo <- tepGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs,lvPlots=FALSE)
	return(tepOutputHandler(res=res,tepPlotInfo=tepPlotInfo))
}
