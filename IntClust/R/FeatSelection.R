FeatSelection<-function(List,Selection=NULL,BinData,ContData=NULL,Datanames=NULL,nrclusters=NULL,topC=NULL,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE){
	
	
	if(is.null(Datanames)){
		for(j in 1:(length(BinData)+length(ContData))){
			Datanames[j]=paste("Data",j,sep=" ")	
		}
	}
	
	if(class(Selection)=="character"){
		
		for(i in 1:length(BinData)){
			BinData[[i]]=BinData[[i]]+0
			BinData[[i]]<-BinData[[i]][,which(colSums(BinData[[i]]) != 0 & colSums(BinData[[i]]) != nrow(BinData[[i]]))]
		}
		
		cpdSet <- rownames(BinData[[1]])
		
		ResultFeat=list()
		Characteristics=list()
		temp=list()
		
		LeadCpds=Selection #names of the compounds
		OrderedCpds=cpdSet
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")
		
		group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
		
		#Determine characteristic features for the compounds: fishers exact test
		result=list()
		for(j in 1: length(BinData)){
			binMat=BinData[[j]]

			pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
			
			pFish <- sort(pFish)
			adjpFish<-stats::p.adjust(pFish, method = "fdr")
			
			AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
			AllFeat$Names=as.character(AllFeat$Names)
			if(is.null(topC)){
				topC=length(which(pFish<sign))
			}
			
			TopFeat=data.frame(Names=as.character(names(pFish[0:topC])),P.Value=pFish[0:topC],adj.P.Val=adjpFish[0:topC])
			TopFeat$Names=as.character(TopFeat$Names)
			temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
			result[[j]]<-temp1
			names(result)[j]=Datanames[j]
			
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
			if(is.null(topC)){
				topC=length(which(pTTest<sign))
			}
			
			TopFeat=data.frame(Names=as.character(names(pTTest[0:topC])),P.Value=pTTest[0:topC],adj.P.Val=adjpTTest[0:topC])
			TopFeat$Names=as.character(TopFeat$Names)
			temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
			resultC[[j]]<-temp1
			names(resultC)[j]=Datanames[length(BinData)+j]
			
		}
		}
		
		temp[[2]]=c(result,resultC)
		
		
		names(temp)=c("Compounds","Characteristics")
		
		ResultFeat[[1]]=temp
		names(ResultFeat)="Selection"
		
		
		
	}
	
	else if(class(Selection)=="numeric" & !(is.null(List))){
		
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
		
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		List=ListNew
		for(i in 1:length(BinData)){
			BinData[[i]]=BinData[[i]]+0
			BinData[[i]]<-BinData[[i]][,which(colSums(BinData[[i]]) != 0 & colSums(BinData[[i]]) != nrow(BinData[[i]]))]
		}
		
		cpdSet <- rownames(BinData[[1]])
		
		ResultFeat=list()
		for(k in 1:dim(Matrix)[1]){	
			cluster=Selection

			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			
			Genes=list()
			temp=list()

			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the compounds
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")
			
			group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
			
			#Determine characteristic features for the compounds: fishers exact test
			result=list()
			for(j in 1: length(BinData)){
				binMat=BinData[[j]]
				pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
				
				pFish <- sort(pFish)
				adjpFish<-stats::p.adjust(pFish, method = "fdr")
				
				AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
				AllFeat$Names=as.character(AllFeat$Names)
				if(is.null(topC)){
					topC=length(which(pFish<0.05))
				}
				
				TopFeat=AllFeat[0:topC,]
				TopFeat$Names=as.character(TopFeat$Names)
				temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
				result[[j]]<-temp1
				names(resultC)[j]=Datanames[length(BinData)+j]
				
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
					if(is.null(topC)){
						topC=length(which(pTTest<sign))
					}
					
					TopFeat=data.frame(Names=as.character(names(pTTest[0:topC])),P.Value=pTTest[0:topC],adj.P.Val=adjpTTest[0:topC])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultC[[j]]<-temp1
					names(resultC)[j]=Datanames[j]
					
				}
			}
			
			temp[[2]]=c(result,resultC)
			
			names(temp)=c("Compounds","Characteristics")
			ResultFeat[[k]]=temp
			names(ResultFeat)[k]=paste(names[k],": Cluster",cluster, sep=" ")
		}		
	}
	
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultFeat)
	
	
	
}