allanVarSelect <-
function(BaseModel,TrnDataSetFile,ValDataSetFile,ResponseCol=1,NumOfSteps=10,criteria="AIC",currentchunksize=-1,silent=TRUE,MemoryAllowed=0.5,TestedRows=1000,AdjFactor=0.095){
#this is the main work function that takes a biglm object and performs a forward variable selection routine on it
#default file to test memory or models is Training data
filename=TrnDataSetFile
#get the optimal chunksize
#get optimal chunksize if not specified
if (currentchunksize==-1){
currentchunksize<-getbestchunksize(filename,MemoryAllowed=MemoryAllowed,TestedRows=TestedRows,AdjFactor=AdjFactor,silent=FALSE)
}

#initialize summary data  
FinalSummary<-NULL
CurrentSummary<-NULL
VarSummary<-NULL

####get all options specified from original or base model###
#to create correct empty model with no vars
BaseCallString <- paste( deparse((BaseModel$call)),collapse="")
OptionStartPos <- gregexpr("," ,BaseCallString)[[1]][1]
OptionString <- substr(BaseCallString,OptionStartPos+1,nchar(BaseCallString))
OptionStartPos2 <- gregexpr("~" ,BaseCallString)[[1]][1]
OptionString2 <- substr(BaseCallString,0,OptionStartPos2)

#get offset option from callstring
OffsetOption <- gregexpr("offset\\(.*?," ,BaseCallString,perl=TRUE)
OffsetOptionStartPos<-OffsetOption[[1]][1]
OffsetOptionEndPos<- (attr(OffsetOption[[1]],"match.length")+OffsetOptionStartPos)
OffsetOptionString <-substr(BaseCallString,OffsetOptionStartPos,OffsetOptionEndPos)

#now look for dataset name and replace with training dataset name
#extracting dataset option name from model
DataOpt1<-gregexpr("data.*?=",OptionString,perl=TRUE)
DataOpt2<-gregexpr("data.*?=.*?,",OptionString,perl=TRUE)
if (DataOpt2[[1]][1]==-1)
{DataOpt2<-gregexpr("data.*?=.*?\\)",OptionString,perl=TRUE)}
DataSetNameStartPos<-DataOpt1[[1]][1]+attr(DataOpt1[[1]],"match.length")[1]
DataSetNameLength<-(attr(DataOpt2[[1]],"match.length")[1]-attr(DataOpt1[[1]],"match.length")[1])
DataSetName<-substr(OptionString,DataSetNameStartPos,DataSetNameStartPos+DataSetNameLength-2)
#assign old data to rawdataset
#OldDataSet<-eval(parse(text=paste(DataSetName)))
#OptionString3<-sub(DataSetName,"datafeed",OptionString,perl=TRUE)

#this is in case an offset is ever used  in glm fitting
#if (OffsetOptionStartPos==-1) 
#{
#NULLModelCall <- parse(text=paste(OptionString2, "NULL ,",OptionString))
#} else {
#NULLModelCall <- parse(text=paste(OptionString2, OffsetOptionString, OptionString))
#}

#create empty model with correct offset and data options
#This is starting point for model
#need to fit a big model on this
columnnames<-names(read.csv(filename, nrows=2,header=TRUE))
ResponseName<-columnnames[ResponseCol]
#first fit smaller model on subsection of data
NULLModelData<-readinbigdata(TrnDataSetFile,chunksize=currentchunksize,col.names=columnnames)
NULLModelData(TRUE)
NULLModelDataSet<-NULLModelData(FALSE)
NULLModelFormula <- parse(text=paste(ResponseName, " ~ NULL"))
OptionStringNULL<-sub(DataSetName,"NULLModelDataSet",OptionString,fixed=TRUE)
NULLModelCall <- parse(text=paste("biglm(", NULLModelFormula ," , ",OptionStringNULL ))
NULLModel<-fitvbiglm(eval(NULLModelCall),TrnDataSetFile,currentchunksize=currentchunksize,silent=silent,MemoryAllowed=MemoryAllowed,TestedRows=TestedRows,AdjFactor=AdjFactor)



#NULLModel<-eval(NULLModelCall)

#create master list of candidate variables to add to model 
VarsToBeAdded<-attr(BaseModel$terms,"term.labels")
NumOrigVars<-length(VarsToBeAdded)

#initialize current model with the null model or intercept model
CurrentModel<-NULLModel

#run the loop that performs variable selection
for (stepCount in 1:NumOfSteps){
#This loop goes through all current variables that still need to be considered
#to enter the model and creates the temporary model
VarMetrics<-NULL
for (m in 1:length(VarsToBeAdded)) {
#create new model with new var
ptext <- parse(text=paste("update.formula(formula(CurrentModel), ~. ","+ ",VarsToBeAdded[m],")"))
TempModelFormulaText <- ptext
#TempModelFormula<-capture.output(eval(ptext))
TempModelFormula<-paste(deparse(eval(ptext)), collapse=" ")
TempModelCall <- parse(text=paste("biglm(", TempModelFormula ," , ",OptionString ))
#fit biglm
TempModel<-fitvbiglm(eval(TempModelCall),TrnDataSetFile,currentchunksize=currentchunksize,silent=silent,MemoryAllowed=MemoryAllowed,TestedRows=TestedRows,AdjFactor=AdjFactor)
#predict on results and get metrics
metricresults<-predictvbiglm(TempModel,ValDataSetFile,currentchunksize=currentchunksize,ResponseCol=ResponseCol,silent=silent,MemoryAllowed=MemoryAllowed,TestedRows=TestedRows,AdjFactor=AdjFactor)
metricresults$Var<-VarsToBeAdded[m]
VarMetrics<-rbind(VarMetrics,metricresults)
}
#save variable selection data for this iteration
VarSummary<-VarMetrics[,1:3]
VarSummary<-as.matrix(VarSummary)
#adds the one which results in minimum MSE  AIC, or  BIC depending on criteria set
if (criteria=="AIC"){
VarToAdd<-which.min(VarSummary[,1])
print("Criteria:AIC")
} else
if (criteria=="BIC"){
VarToAdd<-which.min(VarSummary[,2])
print("Criteria:BIC")
} else
if (criteria=="MSE"){
VarToAdd<-which.min(VarSummary[,3])
print("Criteria:MSE")
} else
{
VarToAdd<-which.min(VarSummary[,1])
print("Criteria flag invalid.  Using AIC")
}

#Change current model.  include last var added
ptext <- parse(text=paste("update.formula(formula(CurrentModel), ~. ","+ ",VarsToBeAdded[VarToAdd],")"))
#UpdatedModelFormula <- capture.output(eval(ptext))
UpdatedModelFormula<-paste(deparse(eval(ptext)), collapse=" ")
UpdatedModelCall <- parse(text=paste("biglm(", UpdatedModelFormula ," , ",OptionString ))
CurrentModel <- eval(UpdatedModelCall)

#delete variable added from list after storing new variable
AddedVar<-VarsToBeAdded[VarToAdd]
VarsToBeAdded<-VarsToBeAdded[-VarToAdd]

#record var added and var statistics
IterationSummary<-(VarSummary[VarToAdd,])
IterationSummary[4]<-AddedVar
IterationSummary[5]<-stepCount

FinalSummary<-rbind(FinalSummary,IterationSummary)
#print old var summary and delete
print("Iteration Result")
print(IterationSummary)
IterationSummary<-NULL
}
CurrentModel$SelectionSummary<-FinalSummary
print("Final Results of Variable Selection:")
print(FinalSummary)
print("Results Stored in $SelectionSummary")
return(CurrentModel)
}

