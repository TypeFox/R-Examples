.distanceheatmaps<-function(Data1,Data2,names,nrclusters=7){
	ClustData1=stats::cutree(Data1,nrclusters) #original clusters (aggregated data clustering)
	ClustData2=stats::cutree(Data2,nrclusters) #clusters of changed method
	
	ClustData1=ClustData1[Data1$order]
	ClustData2=ClustData2[Data2$order]
	
	trueorder1=sort(Data1$order,index.return = TRUE)
	trueorder2=sort(Data2$order,index.return = TRUE)
	
	ordercolors=ClustData1
	order=seq(1,nrclusters)
	
	for (k in 1:length(unique(ClustData1))){
		select=which(ClustData1==unique(ClustData1)[k])
		ordercolors[select]=order[k]
	}
	
	ordercolors2=ClustData2
	
	
	for (k in 1:length(unique(ClustData2))){
		select=which(ClustData2==unique(ClustData2)[k])
		ordercolors2[select]=order[k]
	}
	
	ClustData1=ordercolors[trueorder1$ix]
	ClustData2=ordercolors2[trueorder2$ix]
	
	out=matrix(0,length(ClustData1),length(ClustData2)) #want the rows to be the other method and the columns to be the aggregated data clustering
	
	#names=names[Data2$order]
	
	rownames(out)=names
	colnames(out)=names
	
	for(i in 1:length(names)){
		focus=names[i] #defines the column	
		
		label=ClustData2[i] #color of the cluster is defined by the original cluster that contains focus (1 to 7)		
		
		for(j in 1:length(names)){ #go over the rows
			other=names[j]		
			found=FALSE  #find cluster of other
			k=1
			while(found==FALSE & k<=nrclusters){
				label2=k
				if(other %in% names[ClustData1==label2]){
					othercluster=names[ClustData1==label2]
					found=TRUE
				}
				k=k+1
			}	
			
			if(focus %in% othercluster){ #if other and focus still together: give it color of cluster defined by focus
				out[j,i]=label
			}				
		}				
	}
	return(out)	
}
