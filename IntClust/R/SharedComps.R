SharedComps<-function(List,nrclusters=NULL,fusionsLog=FALSE,WeightClust=FALSE,names=NULL){
	FoundComps=NULL
	
	FoundComps=FindElement("Compounds",List)

	if(is.null(FoundComps)|(is.list(FoundComps) & length(FoundComps) == 0)){
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
		Comps=list()
		for (k in 1:dim(MatrixClusters)[1]){
			clusters=MatrixClusters[k,]

			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			clusters=MatrixClusters[k,]
			Method=list()
			for (i in 1:nrclusters){
			
				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)] 
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds")
				names(temp)="Compounds"
				Method[[i]]=temp
				names(Method)[i]=paste("Cluster",i,sep=" ")
			}
			Comps[[k]]=Method

		}
		names(Comps)=names
		List=Comps
	}
	
	for(a in 1:length(List)){
		if(a==1){
			Comps1=FindElement("LeadCpds",List[[a]])		
		}
		else if(is.list(Comps1) & length(Comps1) != 0){
			Comps2=FindElement("LeadCpds",List[[a]])
			for(b in 1:length(Comps1)){
				Comps1[[b]]=intersect(Comps1[[b]],Comps2[[b]])
			}			
		}
				
	}
	
	namesCl=c()
	for(c in 1:length(Comps1)){
		namesCl=c(namesCl,paste("Cluster",c,sep=" "))
	}
	
	names(Comps1)=namesCl

	
	return(Comps1)
}
