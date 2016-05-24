fitvbiglm <-
function(BaseModel,filename,currentchunksize=-1,silent=TRUE,MemoryAllowed=0.5,TestedRows=1000,AdjFactor=0.095){
#This function opens the connection to the data and fits a linear model.  There is no restriction as
#to whether the dataset can fit into memory

#get optimal chunksize if not specified
if (currentchunksize==-1){
currentchunksize<-getbestchunksize(filename,MemoryAllowed=MemoryAllowed,TestedRows=TestedRows,AdjFactor=AdjFactor,silent=silent)
}

#use the datafeed function to read data
columnnames<-names(read.csv(filename, nrows=2,header=TRUE))
datafeed<-readinbigdata(filename,chunksize=currentchunksize,col.names=columnnames)

#extracting dataset option name from model and replacing it with current
BaseCallString <- paste( deparse((BaseModel$call)),collapse="")
DataOpt1<-gregexpr("data.*?=",BaseCallString,perl=TRUE)
DataOpt2<-gregexpr("data.*?=.*?,",BaseCallString,perl=TRUE)
if (DataOpt2[[1]][1]==-1)
{DataOpt2<-gregexpr("data.*?=.*?\\)",BaseCallString,perl=TRUE)}
DataSetNameStartPos<-DataOpt1[[1]][1]+attr(DataOpt1[[1]],"match.length")[1]
DataSetNameLength<-(attr(DataOpt2[[1]],"match.length")[1]-attr(DataOpt1[[1]],"match.length")[1])
DataSetName<-substr(BaseCallString,DataSetNameStartPos,DataSetNameStartPos+DataSetNameLength-2)
#assign old data to rawdataset
#OldDataSet<-eval(parse(text=paste(DataSetName)))
CallStringNew<-sub(DataSetName,"datafeed(FALSE)",BaseCallString,perl=TRUE)

#fit first interation
#use data grabbing function and initialize for first iteration
datafeed(TRUE)
#call biglm to fit first time
CurrentModel <- eval(parse(text=CallStringNew))
#update rest of the data by using update function
#advance dataset
CurrentDataSet<-datafeed(FALSE)
if (silent!=TRUE){
print("Iterating Through Dataset and Updating Coefficients")
}
while (!is.null(CurrentDataSet)){
CurrentModel<-update(CurrentModel,CurrentDataSet)
CurrentDataSet<-datafeed(FALSE)
}
#return the model
return(CurrentModel)
}

