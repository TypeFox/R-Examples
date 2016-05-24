CompareInteractive<-function(ListM,ListS,nrclusters=NULL,cols=NULL,fusionsLogM=FALSE,fusionsLogS=FALSE,WeightClustM=FALSE,WeightClustS=FALSE,namesM=NULL,namesS=NULL,marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,N=1,...){
	
	MatrixColorsM=ReorderToReference(ListM,nrclusters,fusionsLogM,WeightClustM,namesM)
	
	NamesM=ColorsNames(MatrixColorsM,cols)
	
	nobsM=dim(MatrixColorsM)[2]
	nmethodsM=dim(MatrixColorsM)[1]
	
	if(is.null(namesM)){
		for(j in 1:nmethodsM){
			namesM[j]=paste("Method",j,sep=" ")	
		}
	}
	
	similarM=round(SimilarityMeasure(MatrixColorsM),2)
	nmethodsM=0
	nmethodsS=0
	
	grDevices::dev.new()
	graphics::par(mar=marginsM)
	plotrix::color2D.matplot(MatrixColorsM,cellcolors=NamesM,show.values=FALSE,axes=FALSE,xlab="",ylab="",...)	
	graphics::axis(1,at=seq(0.5,(nobsM-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsM-0.5)),labels=rev(namesM),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsM-0.5)),labels=rev(similarM),cex.axis=0.65,las=2)
	
	if(Interactive==TRUE){
		yseq=c(seq(dim(MatrixColorsM)[1]-0.5,0.5,-1))
		for(i in seq(dim(MatrixColorsM)[1]-0.5,0.5,-1)){
			yseq=c(yseq,rep(i,dim(MatrixColorsM)[2]))
		}
		ids=graphics::identify(x=c(rep(-1,dim(MatrixColorsM)[1]),rep(seq(0.5,dim(MatrixColorsM)[2]-0.5),dim(MatrixColorsM)[1])),y=yseq,n=N,plot=FALSE)
		
		comparison<-function(id){
			if(id%in%seq(dim(MatrixColorsM)[1])){
				grDevices::dev.new()
				graphics::layout(matrix(c(1,2),nrow=2), heights=c(1,2))
				NamesMSel=NamesM[id,]
				namesMSel=namesM[id]
				
				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(t(as.matrix(MatrixColorsM[id,])),cellcolors=NamesMSel,show.values=FALSE,axes=FALSE,xlab="",ylab="")
				graphics::axis(2,at=c(0.5),labels=rev(namesMSel),cex.axis=0.65,las=2)
				#axis(4,at=c(0.5),labels=rev(similarSel),cex.axis=0.65,las=2)
				
				#Find reference for MatrixColorsM
				if(WeightClustM==FALSE){
					temp=FindElement("Results",ListM[1])
					if(!(is.null(temp))&length(temp)!=0){
						Ref=list(Clust=temp$Results_1[[1]])
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else if(length(temp)==0){
						temp=FindElement("Clust",ListM[1])
						Ref=list(Clust=temp$Clust_1)
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else{
						message('Cannot find a reference for the second plot, try: WeightClust=TRUE')
					}
				}
				else{
					Ref=ListM[[1]]
					attr(Ref,'method')<-attributes(ListM[[1]])$method
				}
				
				L=c(Ref,ListS)
				for(i in 1:length(L)){
					if(i==1){
						attr(L[[1]],'method')<-"Ref"
					}
					else{
						attr(L[[i]],"method")<-attributes(ListS[[i-1]])$method
					}					
				}
				MatrixColorsS=ReorderToReference(L,nrclusters,fusionsLogS,WeightClustS,names=c("Ref",namesS))
				MatrixColorsS=MatrixColorsS[-1,]
				NamesS=ColorsNames(MatrixColorsS,cols)
				
				nobs=dim(MatrixColorsS)[2]
				nmethodS=dim(MatrixColorsS)[1]
				
				if(is.null(namesS)){
					for(j in 1:nmethodS){
						namesS[j]=paste("Method",j,sep=" ")	
					}
				}
				
				similarS=round(SimilarityMeasure(MatrixColorsS),2)
				
				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
				graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
				graphics::axis(2,at=c(seq(0.5,nmethodS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
				graphics::axis(4,at=c(seq(0.5,nmethodS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)			
			}
			else{
				SelCluster=t(MatrixColorsM)[id-nrow(MatrixColorsM)]
				Temp=sapply(seq(nrow(MatrixColorsM)),function(i) ncol(MatrixColorsM)*i)
				Row=which(Temp>(id-nrow(MatrixColorsM)))[1]
				Index=which(MatrixColorsM[Row,]!=SelCluster)
				
				grDevices::dev.new()
				graphics::layout(matrix(c(1,2),nrow=2), heights=c(1,2))
				NamesMSel=NamesM[Row,]
				NamesMSel[Index]="white"
				namesMSel=namesM[Row]
				
				graphics::par(mar=marginsM)
				plotrix::color2D.matplot(t(as.matrix(MatrixColorsM[Row,])),cellcolors=NamesMSel,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
				graphics::axis(2,at=c(0.5),labels=rev(namesMSel),cex.axis=0.65,las=2)	
				
				#Find reference for MatrixColorsM
				if(WeightClustM==FALSE){
					temp=FindElement("Results",ListM[1])
					if(!(is.null(temp))&length(temp)!=0){
						Ref=list(Clust=temp$Results_1[[1]])
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else if(length(temp)==0){
						temp=FindElement("Clust",ListM[1])
						Ref=list(Clust=temp$Clust_1)
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else{
						message('Cannot find a reference for the second plot, try: WeightClust=TRUE')
					}
				}
				else{
					Ref=ListM[[1]]
					attr(Ref,'method')<-attributes(ListM[[1]])$method
				}
				
				L=c(Ref,ListS)
				for(i in 1:length(L)){
					if(i==1){
						attr(L[[1]],'method')<-"Ref"
					}
					else{
						attr(L[[i]],"method")<-attributes(ListS[[i-1]])$method
					}					
				}
				
				MatrixColorsS=ReorderToReference(L,nrclusters,fusionsLogS,WeightClustS,names=c("Ref",namesS))
				MatrixColorsS=MatrixColorsS[-1,]
				
				IndexS=lapply(seq(nrow(MatrixColorsS)),function(i) which(MatrixColorsS[i,]!=SelCluster))
				
				NamesS=ColorsNames(MatrixColorsS,cols)
				for(i in 1:nrow(NamesS)){
					NamesS[i,IndexS[[i]]]="white"
				}
				
				nobs=dim(MatrixColorsS)[2]
				nmethods=dim(MatrixColorsS)[1]
				
				if(is.null(namesS)){
					for(j in 1:nmethodsS){
						namesS[j]=paste("Method",j,sep=" ")	
					}
				}
				
				similarS=round(SimilarityMeasure(MatrixColorsS),2)
				
				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
				graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
				graphics::axis(2,at=c(seq(0.5,nmethods-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
				graphics::axis(4,at=c(seq(0.5,nmethods-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)		
				
				
			}
		}
		
		plots=sapply(seq(length(ids)),function(i) comparison(ids[i]))
		
	}
	
}
