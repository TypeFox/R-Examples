epCA <-
function(DATA,DESIGN=NULL,make_design_nominal=TRUE,masses=NULL,weights=NULL,hellinger=FALSE,symmetric=TRUE,graphs=TRUE,k=0){
	main <- deparse(substitute(DATA))
	DATA <- as.matrix(DATA)
	DESIGN<-designCheck(DATA,DESIGN,make_design_nominal)
	
	res <- coreCA(DATA,masses=masses,weights=weights,hellinger=hellinger,symmetric=symmetric,k=k)
	class(res) <- c("epCA","list")
		
	epPlotInfo <- epGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs)
	return(epOutputHandler(res=res,epPlotInfo=epPlotInfo))
}
