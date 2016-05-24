##' @export
`CmptMltPvl.default` <- function (MultiData,CfgFile){
	Conf <- CfgFile
	colnames(Conf) <- c("Contrst","Id","Path") 
	Cntrl <- subset(Conf$Id,Conf$Contrst=="Control")
	Treat <-subset(Conf$Id,Conf$Contrst=="Treat")
	CntrlDt <- MultiData[,1:2]
	CntrlDt$Contrst <- "Cntrl"
	CntrlDt$chr<-NULL
	CntrlDt$start<-NULL
	
	for (ref in Cntrl){
		CntrlDt[[ref]] <- MultiData[[ref]]
	}
	
	CntrlDt<-melt(CntrlDt,na.rm=T, id="Contrst")
	TrTDt <- MultiData[,1:2]
	TrTDt$Contrst <- "Treat"
	TrTDt$chr<-NULL
	TrTDt$start<-NULL
	
	for (trt in Treat){
		TrTDt[[trt]] <- MultiData[[trt]]
	}
	
	TrTDt <- melt(TrTDt,na.rm=T, id="Contrst")
	PrbMlt(rbind(CntrlDt,TrTDt))
}