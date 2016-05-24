ordinalToMetric<-function(data,scaleType="o",patternCoordinates){
if(length(scaleType)==1 && scaleType=="o"){scaleType=rep("o",ncol(data))}
if(length(scaleType)==1 && scaleType=="m"){scaleType=rep("m",ncol(data))}
if(!is.matrix(data) && !is.data.frame(data)){stop("data parameter should be matrix or data frame")}
if(length(scaleType)!=ncol(data)){stop("scale type length should be equal to data column number")}
if(length(scaleType)!=length(patternCoordinates)){stop("scale type length should be equal to data column number")}
for(i in 1:length(scaleType)){
  if(scaleType[i]!="o" & scaleType[i]!="m"){stop("scale type should bo o(ordinal) or m(metric)")}
  if(scaleType[i]=="o" & is.na(patternCoordinates[i])){stop("for scale type o(ordinal) pattern coordinate must be defined")}
}

resul<-NULL;
for(i in 1:ncol(data)){
	if(scaleType[i]=="m"){
		resul<-cbind(resul,data[,i])
	}
	else{
		dVec<-c(data[,i],patternCoordinates[i])
		dim(dVec)<-c(nrow(data)+1,1)
		d<-as.matrix(GDM2(dVec))
		resul<-cbind(resul,1-d[-(nrow(data)+1),(nrow(data)+1)])
	}
}
return (resul)
}