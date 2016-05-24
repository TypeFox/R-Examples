calEdgeCorScore <- function(dataset, class.labels,controlcharactor, edgesbackgrand) {
#This function calculate differential Mutual information
#Inputs:
#     dataset: Marix of gene expression values (rownames are genes, columnnames are samples) 
#     class.labels: Vector of binary labels.
#     controlcharactor: Charactor of control in the class labels.
#     edgesbackgrand: Marix of the edges' backgrand   
#Outputs:
#     EdgeCorScore:Vector of the aberrant correlation in phenotype P based on mutual information (MI) for each edge
#        
       location<-matrix(0,length(edgesbackgrand[,1]),2)
	   location[,1]<-match(edgesbackgrand[,1],rownames(dataset))
	   location[,2]<-match(edgesbackgrand[,2],rownames(dataset))
	   location<-na.omit(location)
	   controlloca<-which(class.labels==controlcharactor)
	   dataset.1<-dataset[location[,1],]
	   dataset.2<-dataset[location[,2],]
	   Cexpress.1<-dataset[location[,1],controlloca]
	   Cexpress.2<-dataset[location[,2],controlloca]
	   EdgeCorScore<-c()
	   for(i in 1:length(location[,1])){
		   EdgeCorScore[i]<-knnmi(dataset.1[i,],dataset.2[i,])-knnmi(Cexpress.1[i,],Cexpress.2[i,])
		  }   
		EdgeID<-matrix(0,length(location[,1]),2)
		EdgeID[,1]<-row.names(dataset.1)
		EdgeID[,2]<-row.names(dataset.2)
		colnames(EdgeID)<-c("Edge1","Edge2")
		Paste<-function(x,c1,c2) paste(x[c1],x[c2],sep="|")
		names(EdgeCorScore)<-apply(EdgeID,1,Paste,c1="Edge1",c2="Edge2")
		return(EdgeCorScore)
	   }