#tepBADA <- function(DATA,scale=TRUE,center=TRUE,DESIGN=NULL,make_design_nominal=TRUE,group.masses=NULL,ind.masses=NULL,weights=NULL,graphs=TRUE,k=0){
tepBADA <- function(DATA,scale=TRUE,center=TRUE,DESIGN=NULL,make_design_nominal=TRUE,group.masses=NULL,weights=NULL,graphs=TRUE,k=0){	
		
	
	DESIGN <- texpoDesignCheck(DATA,DESIGN,make_design_nominal,force_bary=TRUE)
	colDESIGN <- colnames(DESIGN)
	#massedDESIGN<-t(t(DESIGN) * (1/(colSums(DESIGN))))
	massedDESIGN <- t(apply(DESIGN,1,'/',colSums(DESIGN)))
	colnames(massedDESIGN) <- colDESIGN	
	
	main <- deparse(substitute(DATA))		
	DATA <- as.matrix(DATA)

	#XMW <- computeMW(DATA,masses=ind.masses,weights=weights)

	R <- expo.scale(t(massedDESIGN) %*% DATA,scale=scale,center=center)
	this.center <- attributes(R)$`scaled:center`
	this.scale <- attributes(R)$`scaled:scale`	
	RMW <- computeMW(R,masses=group.masses,weights=weights)

	colnames(R) <- colnames(DATA)
	rownames(R) <- colnames(DESIGN)	
	Rdesign <- diag(nrow(R))
	rownames(Rdesign) <- rownames(R)	
	
	#res <- corePCA(R,M=RMW$M,W=RMW$W,k=k)
	res <- epGPCA(R, DESIGN=Rdesign, make_design_nominal=FALSE, scale = FALSE, center = FALSE, masses = RMW$M, weights = RMW$W, graphs = FALSE, k = k)
	res <- res$ExPosition.Data
	res$center <- this.center
	res$scale <- this.scale

	supplementaryRes <- supplementaryRows(DATA,res)
	res$fii <- supplementaryRes$fii
	res$dii <- supplementaryRes$dii
	res$rii <- supplementaryRes$rii
	
	res$lx <- res$fii
	res$ly <- supplementaryCols(t(massedDESIGN),res,center=FALSE,scale=FALSE)$fjj

	assignments <- fii2fi(DESIGN,res$fii,res$fi)
	#assignments$r2 <- R2(RMW$M,res$di,XMW$M,res$dii)
	assignments$r2 <- R2(RMW$M,res$di,ind.masses=NULL,res$dii)
	class(assignments) <- c("tepAssign","list")
	res$assign <- assignments

	#new res here
	class(res) <- c("tepBADA","list")		
	tepPlotInfo <- tepGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs,lvPlots=FALSE)
	return(tepOutputHandler(res=res,tepPlotInfo=tepPlotInfo))
}
