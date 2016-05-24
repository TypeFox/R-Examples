epPCA <-
function(DATA,scale=TRUE,center=TRUE,DESIGN=NULL,make_design_nominal=TRUE,graphs=TRUE,k=0){

	main <- deparse(substitute(DATA))
	DESIGN <- designCheck(DATA,DESIGN,make_design_nominal)
	DATA <- as.matrix(DATA)
	DATA <- expo.scale(DATA,scale=scale,center=center)
	this.center <- attributes(DATA)$`scaled:center`
	this.scale <- attributes(DATA)$`scaled:scale`

	res <- corePCA(DATA,k=k)
	res$center <- this.center
	res$scale <- this.scale
	class(res) <- c("epPCA","list")
	
	epPlotInfo <- epGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs)
	return(epOutputHandler(res=res,epPlotInfo=epPlotInfo))
}
