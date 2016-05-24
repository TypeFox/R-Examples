CharacteristicFeatures<-function(List,Selection=NULL,BinData,ContData=NULL,Datanames=NULL,nrclusters=NULL,sign=0.05,topC=NULL,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	if(is.null(Datanames)){
		for(j in 1:(length(BinData)+length(ContData))){
			Datanames[j]=paste("Data",j,sep=" ")	
		}
	}
	
	if(!(is.null(Selection))){
		ResultFeat=FeatSelection(List,Selection,BinData,ContData,Datanames,nrclusters,topC,sign,fusionsLog,WeightClust)
	}
	else{
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(WeightClust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		
		List=ListNew
		for(i in 1:length(BinData)){
			BinData[[i]]=BinData[[i]]+0
			BinData[[i]]<-BinData[[i]][,which(colSums(BinData[[i]]) != 0 & colSums(BinData[[i]]) != nrow(BinData[[i]]))]
		}
		
		cpdSet <- rownames(BinData[[1]])
		
		ResultFeat=list()
		maxclus=0
		for (k in 1:dim(MatrixClusters)[1]){

			clusters=MatrixClusters[k,]
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			Characteristics=list()
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			for (i in clust){		
				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)] 
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds") #names of the compounds
				
				group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
				
				#Determine characteristic features for the compounds: fishers exact test
				result=list()
				for(j in 1: length(BinData)){
					binMat=BinData[[j]]
									
					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")
					
					AllFeat=data.frame(Names=names(pFish),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)
					
					if(is.null(topChar)){
						topChar=length(which(pFish<sign))
					}
					
					TopFeat=AllFeat[0:topChar,]
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(result)[j]=Datanames[length(BinData)+j]
					
				}
				
				
				resultC=list()
				if(!is.null(ContData)){
					for(j in 1:length(ContData)){
						contMat=ContData[[j]]
						
						group1=which(group==1)
						group2=which(group==0)
						
						
						pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)
						
						pTTest <- sort(pTTest)
						adjpTTest<-stats::p.adjust(pTTest, method = "fdr")
						
						AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
						AllFeat$Names=as.character(AllFeat$Names)
						if(is.null(topChar)){
							topChar=length(which(pTTest<sign))
						}
						
						TopFeat=data.frame(Names=as.character(names(pTTest[0:topChar])),P.Value=pTTest[0:topChar],adj.P.Val=adjpTTest[0:topChar])
						TopFeat$Names=as.character(TopFeat$Names)
						temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
						resultC[[j]]<-temp1
						names(resultC)[j]=Datanames[j]
						
					}
				}
				
				temp[[2]]=c(result,resultC)
				
				names(temp)=c("Compounds","Characteristics")
				
				Characteristics[[i]]=temp
				
				names(Characteristics)[i]=paste("Cluster",i,sep=" ")
			}
			ResultFeat[[k]]=Characteristics
			
		}
		names(ResultFeat)=names
		for(i in 1:length(ResultFeat)){
			for(k in 1:length(ResultFeat[[i]])){
				if(is.null(ResultFeat[[i]][[k]])[1]){
					ResultFeat[[i]][[k]]=NA
					names(ResultFeat[[i]])[k]=paste("Cluster",k,sep=" ")
				}			
			}
			if(length(ResultFeat[[i]]) != maxclus){
				extra=maxclus-length(ResultFeat[[i]])
				for(j in 1:extra){
					ResultFeat[[i]][[length(ResultFeat[[i]])+1]]=NA
					names(ResultFeat[[i]])[length(ResultFeat[[i]])]=paste("Cluster",length(ResultFeat[[i]]),sep=" ")
				}
			}
		} 	
		
	}
		
	return(ResultFeat)	
}
