CompareSvsM<-function(ListS,ListM,nrclusters=NULL,cols=NULL,fusionsLogS=FALSE,fusionsLogM=FALSE,WeightClustS=FALSE,WeightClustM=FALSE,namesS=NULL,namesM=NULL,margins=c(8.1,3.1,3.1,4.1),plottype="new",location=NULL,...){
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new(wdith=14,height=7)
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	nmethodsS=0
	nmethodsM=0
	
	MatrixColorsS=ReorderToReference(ListS,nrclusters,fusionsLogS,WeightClustS,namesS)
	MatrixColorsM=ReorderToReference(c(ListS[1],ListM),nrclusters,fusionsLogM,WeightClustM,c("ref",namesM))
	
	similarS=round(SimilarityMeasure(MatrixColorsS),2)
	similarM=round(SimilarityMeasure(MatrixColorsM),2)
	
	MatrixColorsM=MatrixColorsM[-c(1),]
	
	NamesM=ColorsNames(MatrixColorsM,cols)
	NamesS=ColorsNames(MatrixColorsS,cols)
	
	nobsM=dim(MatrixColorsM)[2]
	nmethodsM=dim(MatrixColorsM)[1]
	
	nobsS=dim(MatrixColorsS)[2]
	nmethodsS=dim(MatrixColorsS)[1]
	
	if(is.null(namesS)){
		for(j in 1:nmethodsS){
			namesS[j]=paste("Method",j,sep=" ")	
		}
	}
	
	if(is.null(namesM)){
		for(j in 1:nmethodsM){
			namesM[j]=paste("Method",j,sep=" ")	
		}
	}
		

	plottypein(plottype,location)
	graphics::par(mfrow=c(1,2),mar=margins)
	plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="",...)
	graphics::axis(1,at=seq(0.5,(nobsS-0.5)),labels=colnames(MatrixColorsS),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)
	
	plotrix::color2D.matplot(MatrixColorsM,cellcolors=NamesM,show.values=FALSE,axes=FALSE,xlab="",ylab="",...)
	graphics::axis(1,at=seq(0.5,(nobsM-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsM-0.5)),labels=rev(namesM),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsM-0.5)),labels=rev(similarM[-1]),cex.axis=0.65,las=2)
	plottypeout(plottype)

}


