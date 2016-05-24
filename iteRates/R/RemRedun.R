RemRedun <-
function(ClusRes){
	for(G in 1:length(ClusRes)){
	i=1
	 tempC<-ClusRes[[G]]$group
	 tempV<-ClusRes[[G]]$VecBin
	 tempM<-array(dim=c(length(tempC[,1]),length(tempV)))
	 	for (j in 1:length(tempC[,1])){
	 		for (k in 1:length(tempV)){
	 			grpend<-cumsum(tempV)
	 			grestart<-c(1,1+grpend[-length(grpend)])
	 				pasted<-paste(as.character(tempV[k]),paste(tempC[j,grestart[k]:grpend[k]],collapse=""),sep="")
	 				tempM[j,k]<-pasted
	 			}
	 		
	 		}
		tempM2<-t(apply(tempM,1,sort))
		pasted<-NA
			for (i in 1:length(tempM2[,1])){
				pasted[i]<-paste(tempM2[i,],collapse="")
				}
		uclus<-unique(pasted)
			for (i in 1:length(uclus)){
				temp2exclue<-which(uclus[i]==pasted)
				if(length(temp2exclue)>1){
					for(j in 2:length(temp2exclue)){
						tempC[temp2exclue[j],]<-NA
				
						
					}
				}
			tempCna<-as.matrix(na.exclude(as.data.frame(tempC)))
		
			ClusRes[[G]]$group<-as.array(tempCna)
			ClusRes[[G]]$dimen<-dim(tempCna)
				}
			}
		return(ClusRes)	
		}

