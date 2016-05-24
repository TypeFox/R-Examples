tepGPLS <- function(DATA1,DATA2,center1=TRUE,scale1="SS1",center2=TRUE,scale2="SS1",DESIGN=NULL,make_design_nominal=TRUE,weights1=NULL,weights2=NULL,graphs=TRUE,k=0){

	if(nrow(DATA1) != nrow(DATA2)){
		stop("DATA1 and DATA2 must have the same number of rows.")
	}

	main <- paste("GPLS: ",deparse(substitute(DATA1))," & ", deparse(substitute(DATA2)),sep="")
	DESIGN <- texpoDesignCheck(DATA1,DESIGN,make_design_nominal=make_design_nominal)
	DESIGN <- texpoDesignCheck(DATA2,DESIGN, make_design_nominal=FALSE)	

	DATA1 <- as.matrix(DATA1)
	DATA2 <- as.matrix(DATA2)	
	DATA1 <- expo.scale(DATA1,scale=scale1,center=center1)	
	DATA2 <- expo.scale(DATA2,scale=scale2,center=center2)		
	MW1 <- computeMW(DATA1,weights=weights1)
	MW2 <- computeMW(DATA2,weights=weights2)

	R <- t(DATA1) %*% DATA2

	#res <- corePCA(R,M=MW1$W,W=MW2$W,k=k)
	res <- epGPCA(DATA=R,masses=MW1$W,weights=MW2$W,k=k,graphs=FALSE,scale=FALSE,center=FALSE)
	res <- res$ExPosition.Data
	res$center <- NULL
	res$scale <- NULL	
	res$W1 <- res$M
	res$W2 <- res$W
	res$M <- res$W <- NULL
	res$data1.norm <- list(center=attributes(DATA1)$`scaled:center`,scale=attributes(DATA1)$`scaled:scale`)
	res$data2.norm <- list(center=attributes(DATA2)$`scaled:center`,scale=attributes(DATA2)$`scaled:scale`)	
	
	##we should probably include the weights into the factor scores across all methods.
	res$lx <- supplementalProjection(DATA1,res$fi * matrix(res$W1,nrow(res$fi),ncol(res$fi),byrow=TRUE),Dv=res$pdq$Dv)$f.out
	res$ly <- supplementalProjection(DATA2,res$fj * matrix(res$W2,nrow(res$fj),ncol(res$fj),byrow=TRUE),Dv=res$pdq$Dv)$f.out

	class(res) <- c("tepGPLS","list")	
	tepPlotInfo <- tepGraphs(res=res,DESIGN=DESIGN,main=main,graphs=graphs)	
	return(tepOutputHandler(res=res,tepPlotInfo=tepPlotInfo))
}