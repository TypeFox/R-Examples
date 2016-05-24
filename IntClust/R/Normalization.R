Normalization<-function(Data,method=c("Quantile","Fisher-Yates","Standardize","Range","Q","q","F","f","S","s","R","r")){
		
		method=substring(method[1],1,1)  #Function also accepts begin letter of each method	
		method=match.arg(method)
		if(method=="S"|method=="s"){
			Data1<-pls::stdize(Data)#center and divide by standard deviation
			return(Data1)
		}
		
		else if(method=="R"|method=="r"){
#			rangenorm<-function(x){
#				minc=min(x)
#				maxc=max(x)
#				rx=(x-minc)/(maxc-minc)
#			}
#			
#			DataN=apply(Data,2,rangenorm)
#			Data=DataN
#			return(Data)
		
			minc=min(as.vector(Data))
			maxc=max(as.vector(Data))
			DataN=(Data-minc)/(maxc-minc)
			Data=DataN
			return(Data)

		}
		else if(method=="Q"|method=="q"){
			DataColSorted <- apply(Data, 2, sort)
			NValues=apply(DataColSorted,1,stats::median,na.rm=T)
		}
		else{
			DataColSorted <- apply(Data, 2, sort)
			NValues=stats::qnorm(1:nrow(Data)/(nrow(Data)+1))
		}
		
	DataRanked=c(apply(Data,2,rank))
		
	array(stats::approx(1:nrow(Data),NValues,DataRanked)$y,dim(Data),dimnames(Data))
	t=stats::approx(1:nrow(Data),NValues,DataRanked)
	return(Data)
	
}

