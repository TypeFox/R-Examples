##' @export

`ProcMutilDt.default` <- function (MultiData){
	colnames(MultiData$Conf) <- c("Contrst","Id","Path") 
	Cntrl <- subset(MultiData$Conf$Id,MultiData$Conf$Contrst=="Control")
	Treat <-subset(MultiData$Conf$Id,MultiData$Conf$Contrst=="Treat")
	TmpCntrl <- MultiData$MultiSvScl[,1:3]
	
	if(length(Cntrl) >1){
		for (ref in Cntrl){
			TmpCntrl[[ref]] <- MultiData$MultiSvScl[[ref]]
		}
		TmpCntrl$CntrlMean <- rowMeans(TmpCntrl[,4:dim(TmpCntrl)[2]])
	}else{
  		for (ref in Cntrl){
    		TmpCntrl[[ref]] <- MultiData$MultiSvScl[[ref]]
  		}
  	TmpCntrl$CntrlMean <- TmpCntrl[,4]}
	TmpMultiSVOutLog2 <- MultiData$MultiSvScl[,1:3]
	
	for (sample in Treat){
		TmpMultiSVOutLog2[[sample]] <- log2(MultiData$MultiSvScl[[sample]]/TmpCntrl$CntrlMean) 
	}
	TmpMultiSVOutLog2$PopLg2 <- rowMeans(TmpMultiSVOutLog2[4:dim(TmpMultiSVOutLog2)[2]])
  MultiData$MultiSvLg2 <- TmpMultiSVOutLog2;rm(TmpMultiSVOutLog2)
  MultiData
}

