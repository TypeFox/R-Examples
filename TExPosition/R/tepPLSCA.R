#Hellinger not available yet. Not until I can get the caNorm pipeline fixed.
tepPLSCA <-
function(DATA1,DATA2,make_data1_nominal=FALSE,make_data2_nominal=FALSE,DESIGN=NULL,make_design_nominal=TRUE,weights1=NULL,weights2=NULL,symmetric=TRUE,graphs=TRUE,k=0){

	main <- paste("PLSCA: ",deparse(substitute(DATA1))," & ", deparse(substitute(DATA2)),sep="")		
	if(nrow(DATA1) != nrow(DATA2)){
		stop("DATA1 and DATA2 must have the same number of rows.")
	}
	if(make_data1_nominal){
		DATA1 <- makeNominalData(DATA1)
	}
	if(make_data2_nominal){
		DATA2 <- makeNominalData(DATA2)
	}	
	
	DESIGN <- texpoDesignCheck(DATA1,DESIGN,make_design_nominal)
	DESIGN <- texpoDesignCheck(DATA2,DESIGN,make_design_nominal=FALSE)	

	DATA1 <- as.matrix(DATA1)
	DATA2 <- as.matrix(DATA2)	
	R <- t(DATA1) %*% DATA2
	
	res <- coreCA(R,masses=weights1,weights=weights2,hellinger=FALSE,symmetric=symmetric,k=k)
	#res <- epCA(R, masses = weights1, weights = weights2, hellinger = FALSE, symmetric = symmetric, graphs = FALSE, k=k)
	#res <- res$ExPosition.Data
	res$W1 <- res$M
	res$W2 <- res$W
	res$M <- res$W <- NULL
	res$lx <- supplementalProjection(makeRowProfiles(DATA1)$rowProfiles,res$fi,Dv=res$pdq$Dv)$f.out / sqrt(nrow(DATA1))
	if(symmetric){
		res$ly <- supplementalProjection(makeRowProfiles(DATA2)$rowProfiles,res$fj,Dv=res$pdq$Dv)$f.out / sqrt(nrow(DATA2))
	}else{
		res$ly <- (supplementalProjection(makeRowProfiles(DATA2)$rowProfiles,res$fj,Dv=rep(1,length(res$pdq$Dv)))$f.out / sqrt(nrow(DATA2)))
	}
	
	class(res) <- c("tepPLSCA","list")
	tepPlotInfo <- tepGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs)	
	return(tepOutputHandler(res=res,tepPlotInfo=tepPlotInfo))
}