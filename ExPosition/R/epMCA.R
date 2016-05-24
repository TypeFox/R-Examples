epMCA <-
function(DATA,make_data_nominal=TRUE,DESIGN=NULL,make_design_nominal=TRUE,masses=NULL,weights=NULL,hellinger=FALSE,symmetric=TRUE,correction=c('b'),graphs=TRUE,k=0){
	main <- deparse(substitute(DATA))	
	DATA <- as.matrix(DATA)
	if(make_data_nominal){
		nominalData <- makeNominalData(DATA)
	}else{
		nominalData <- nominalCheck(DATA)
	}
	
	DESIGN<-designCheck(DATA,DESIGN,make_design_nominal)

	res <- coreCA(nominalData,masses=masses,weights=weights,hellinger=hellinger,symmetric=symmetric,k=k)
	
	##if they give the wrong correction
	weCanCorrect <- ((('b' %in% correction) || ('g' %in% correction)) && !hellinger)
	##if the data table is nominal, the col & row Sums should be identical and have 1 unique element.
	numUniqueRows <- length(unique(round(rowSums(nominalData))))
	if( weCanCorrect && (numUniqueRows==1) ){
		res <- mca.eigen.fix(nominalData,res, make_data_nominal=FALSE,correction=correction,symmetric=symmetric)
	}else{
		res$pdq.uncor <- res$pdq
	}
	class(res) <- c("epMCA","list")	
	
	epPlotInfo <- epGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs)
	return(epOutputHandler(res=res,epPlotInfo=epPlotInfo))
}
