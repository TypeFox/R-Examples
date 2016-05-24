rankInverseNormalDataFrame <-
function(variableList,data,referenceframe) 
{
  
	colnamesList <- as.vector(variableList[,1]);
	size = length(colnamesList);
	SortedCtr = referenceframe;
	nrowsCtr =  nrow(referenceframe);
	nrows = nrow(data);
	zRankInverseFrame = data;
	for (i in 1:size)
	{ 
		if (length(table(data[,colnamesList[i]]))>2)
		{
			SortedCtr[,colnamesList[i]]<-SortedCtr[order(referenceframe[,colnamesList[i]]),colnamesList[i]];
			minvalue <- SortedCtr[1,colnamesList[i]];
			maxvalue <- SortedCtr[nrowsCtr,colnamesList[i]];
			cat(" Variable: ",colnamesList[i],"Min: ",minvalue," Max: ",maxvalue,"\n")      
			zRankInverseFrame[,colnamesList[i]]<-.Call("rankInverseNormalCpp",nrows,data[,colnamesList[i]],minvalue,maxvalue,SortedCtr[,colnamesList[i]]);
      
		}
	}
	return (zRankInverseFrame);
}
