predictvbiglm <-
function(BaseModel,ValFileName,currentchunksize=-1,ResponseCol=1,silent=TRUE,MemoryAllowed=0.5,TestedRows=1000,AdjFactor=0.095){
#This function takes a biglm object and makes predictions as well as returns fit statistics 

#get optimal chunksize if not specified
if (currentchunksize==-1){
currentchunksize<-getbestchunksize(ValFileName,MemoryAllowed=MemoryAllowed,TestedRows=TestedRows,AdjFactor=AdjFactor,silent=silent)

}

#this function fills empty rows with NAs
fillrows<-function(InData,fulllength){
#InData<-RawPredictions
#fulllength<-CurrentNumObs
rowmat<-data.matrix(rownames(InData))
#rownames(rowmat)<-rownames(InData)
allrows<-matrix(1:fulllength)
matched<-data.matrix(match(allrows,rowmat))
rownames(matched)<-allrows
unmatched<-subset(matched,is.na(matched)==TRUE)
#now combine with original and sort by rowname
TempData<-rbind(InData,unmatched)
OutData<-TempData[as.character(sort(as.numeric(rownames(TempData)))),]
OutData<-as.matrix(OutData)
return(OutData)
}

filterweights<-function(InWeights,InData,fulllength){
#InWeights<-weightvector
#InData<-RawPredictions
#fulllength<-RawCurrentNumObs
rowmat<-data.matrix(rownames(InData))
#rownames(rowmat)<-rownames(InData)
allrows<-matrix(1:fulllength)
matched<-data.matrix(match(allrows,rowmat))
rownames(matched)<-allrows
unmatched<-subset(matched,is.na(matched)==TRUE)
#now keep only valid weightvect values
OutWeights<-data.matrix(InWeights[as.numeric(rowmat)])
rownames(OutWeights)<-rownames(InData)
#now combine with original and sort by rowname
TempData<-rbind(OutWeights,unmatched)
OutData<-TempData[as.character(sort(as.numeric(rownames(TempData)))),]
OutData<-as.matrix(OutData)
return(OutData)

}



#get datafeed
columnnames<-names(read.csv(ValFileName, nrows=2,header=TRUE))
datafeed<-readinbigdata(ValFileName,chunksize=currentchunksize,col.names=columnnames)

#initialize running total variables
ObsVec<-NULL
RSSVec<-NULL
VarianceVec<-NULL
MeanVec<-NULL
Var1Vec<-NULL
YValues<-NULL
WeightVec<-NULL

#fit first interation
#use data grabbing function and initialize for first iteration
datafeed(TRUE)
#use predict for biglm to get predictions
CurrentDataSet<-datafeed(FALSE)
while (!is.null(CurrentDataSet)){
RawPredictions<-predict(BaseModel,CurrentDataSet)

#set current number of observations in current iteration
RawCurrentNumObs<-dim(CurrentDataSet)[1]
NewCurrentNumObs<-dim(RawPredictions)[1]



#now fill in predictions that were not predicted with NAs
#check first to see if neceesary
if (RawCurrentNumObs==NewCurrentNumObs){
Predictions<-RawPredictions
CurrentNumObs<-dim(Predictions)[1]
}
else{
Predictions<-fillrows(RawPredictions,RawCurrentNumObs)
CurrentNumObs<-dim(Predictions)[1]
}


#assign weight vector depending on whether weights have been specified
#weights not assigned
if (is.null(BaseModel$weights)){
#assign weights as all same
weightvector<-as.vector(matrix(1,CurrentNumObs,1))
} 
#weights assigned
else{
#parse name of weights
weightname<-substr(BaseModel$weights,1,10000)[2]
#assign to weights
tempweightvector<-as.vector(eval(parse(text=paste("CurrentDataSet","$",weightname,sep=""))))
#now filter out weights with obs that were not used if needed
if (RawCurrentNumObs==NewCurrentNumObs){
weightvector<-filterweights(tempweightvector,RawPredictions,RawCurrentNumObs)
 }else{
weightvector<-tempweightvector
}
}

#calculate mean, variance, and RSS for current chunk
CurrentMean=weighted.mean(CurrentDataSet[,ResponseCol],weightvector,na.rm=TRUE)
CurrentVariance=sum(((CurrentDataSet[,ResponseCol]-CurrentMean)^2)*weightvector,na.rm=TRUE)/sum(weightvector,na.rm=TRUE)
CurrentRSS= sum(((Predictions-CurrentDataSet[,ResponseCol])^2)*weightvector,na.rm=TRUE)

#store data for all chunks in vector for final calculation
RSSVec<-c(RSSVec,CurrentRSS)
MeanVec<-c(MeanVec,CurrentMean)
VarianceVec<-c(VarianceVec,CurrentVariance)
ObsVec<-c(ObsVec,CurrentNumObs)
WeightVec<-c(WeightVec,weightvector)
#for debugging save y values
#YValues<-c(YValues,CurrentDataSet[,ResponseCol])

#increment dataset feed connection to next chunk
CurrentDataSet<-datafeed(FALSE)
}

#calculate total RSS mean and variance
TotalRSS<-sum(RSSVec,na.rm=TRUE)
TotalMean<-weighted.mean(MeanVec,ObsVec,na.rm=TRUE)
TotalObs<-sum(ObsVec,na.rm=TRUE)
MSEValue<-TotalRSS/sum(WeightVec,na.rm=TRUE)

MeanVariance<-sum((MeanVec)^2*(ObsVec/TotalObs),na.rm=TRUE)-TotalMean^2
TotalVariance<-sum((VarianceVec*ObsVec)/TotalObs,na.rm=TRUE)+MeanVariance
#Calculate degrees of freedom: add all variables plus intercept if necessary
NumParameters=attr(BaseModel$terms,"intercept")+length(attr(BaseModel$terms,"term.labels"))

#Calculate likelihood and resulting measures
LogLikelihood= (-TotalObs/2)*(log(2*pi*(TotalVariance^2)))-(1/2)*TotalRSS/(TotalVariance^2)
AICValue=2*NumParameters-2*LogLikelihood
BICValue=-2*LogLikelihood+NumParameters*log(TotalObs)
#data frame with values
metricvalues<-data.frame(AICValue,BICValue,MSEValue)
#return metric for measuring model
return(metricvalues)

}

