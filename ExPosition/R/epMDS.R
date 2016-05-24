epMDS <-
function(DATA,DATA_is_dist=TRUE,method="euclidean",DESIGN=NULL,make_design_nominal=TRUE,masses=NULL,graphs=TRUE,k=0){
	main <- deparse(substitute(DATA))	
	DESIGN<-designCheck(DATA,DESIGN,make_design_nominal)
	DATA <- as.matrix(DATA)
		
	if(DATA_is_dist && (nrow(DATA)==ncol(DATA))){
		D <- DATA
		MW <- computeMW(D,masses=masses)
	}else{
		#print('Creating distance matrix from DATA.')
		D.MW <- makeDistancesAndWeights(DATA,method=method,masses=masses)
		D <- D.MW$D #already squared from the above function.
		MW <- D.MW$MW		
	}
	#S <- mdsTransform(D,MW)
	#res <- coreMDS(S,MW$M,k=k)
	
	res <- coreMDS(D,MW$M,k=k)
	#overwrite res with half of res because it's MDS and we don't care
	res <- list(fi=res$fi,di=res$di,ci=res$ci,ri=res$ri,t=res$t,eigs=res$eigs,pdq=res$pdq,M=res$masses,X=res$X,D=D)
	class(res) <- c("epMDS","list")

	epPlotInfo <- epGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs)
	return(epOutputHandler(res=res,epPlotInfo=epPlotInfo))
}
