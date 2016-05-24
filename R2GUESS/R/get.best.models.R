get.best.models <- function(path.out,path.in=NULL,file.in,label.X=NULL,p,MAP.file=NULL){

#p number of VARIABLE
NameBestModel <- file.path(path.out, paste(file.in,"output_best_visited_models.txt",sep="_"))

if(is.null(label.X)) label.X <- 1:p

#Creating the Best Model object
if(is.null(MAP.file)){
BestModels<-list('Rank'=c(),'nVisits'=c(),'FirstVisit'=c(),'nEvalBefore1st'=c(),'ModeSize'=c(),'logCondPost'=c(),'postProb'=c(),'jeffries'=c(),'modelName'=c())}else{
if(is.data.frame(MAP.file)){VARIABLELabels <- MAP.file}else{
NameMap.file <- file.path(path.in, MAP.file)
VARIABLELabels <- read.table(NameMap.file,header=TRUE)}

BestModels<-list('Rank'=c(),'nVisits'=c(),'FirstVisit'=c(),'nEvalBefore1st'=c(),'ModeSize'=c(),'logCondPost'=c(),'postProb'=c(),'jeffries'=c(),'modelPosInX'=c(),'modelCHR'=c(),'modelPosn'=c(),'modelName'=c())

}


#Reading the file
TmpRead <- readLines(NameBestModel)
#Removing first line: headers
TmpRead <- TmpRead[-1]
MyLine <- 1
Stop <- 0
while(Stop==0){
  tmpRow<-unlist(strsplit(TmpRead[MyLine],"\t"))
  BestModels$Rank[MyLine] <- as.numeric(tmpRow[1])
  BestModels$nVisits[MyLine] <- as.numeric(tmpRow[2])
  BestModels$FirstVisit[MyLine] <- as.numeric(tmpRow[3])
  BestModels$nEvalBefore1st[MyLine] <- as.numeric(tmpRow[4])
  nVARIABLEin <- as.numeric(tmpRow[5])
  BestModels$ModeSize[MyLine] <- nVARIABLEin
  BestModels$logCondPost[MyLine] <- as.numeric(tmpRow[6])
  BestModels$postProb[MyLine] <- as.numeric(tmpRow[7])
  BestModels$jeffries[MyLine] <- as.numeric(tmpRow[8])
		
if(nVARIABLEin>0){
    for(VARIABLE in 1:nVARIABLEin){
      tmpCurrVARIABLE <- as.numeric(tmpRow[8+VARIABLE])
if(!is.null(MAP.file)){
    tmpCurrCHR <- as.numeric(VARIABLELabels$Chr[tmpCurrVARIABLE])
      tmpCurrPosn <- as.numeric(VARIABLELabels$Posn[tmpCurrVARIABLE])
	tmpCurrName <- as.character(VARIABLELabels$SNPName[tmpCurrVARIABLE])}else{
      tmpCurrName <- as.character(label.X[tmpCurrVARIABLE])
}

      if(VARIABLE==1){
	if(!is.null(MAP.file)){
       CurrVARIABLE <- tmpCurrVARIABLE
        CurrCHR <- tmpCurrCHR
        CurrPosn <- tmpCurrPosn}
        CurrName <- tmpCurrName
      }
      if(VARIABLE>1){
	if(!is.null(MAP.file)){
        CurrVARIABLE <- paste(CurrVARIABLE,tmpCurrVARIABLE,sep=' ')
        CurrCHR <- paste(CurrCHR,tmpCurrCHR,sep=' ')
        CurrPosn <- paste(CurrPosn,tmpCurrPosn,sep=' ')}
        CurrName <- paste(CurrName,tmpCurrName,sep=' ')

      }
    }
  }
  if(nVARIABLEin==0){
if(!is.null(MAP.file)){
    CurrVARIABLE=NA
    CurrCHR=NA
    CurrPosn=NA}
    CurrName=NA
  }
if(!is.null(MAP.file)){
  BestModels$modelPosInX[MyLine]=CurrVARIABLE
  BestModels$modelCHR[MyLine]=CurrCHR
  BestModels$modelPosn[MyLine]=CurrPosn}
  BestModels$modelName[MyLine]=CurrName
#  print(MyLine)
  MyLine <- MyLine+1

  if(MyLine>length(TmpRead)){
    Stop <- 1
  }
}

ModelsWriteName <- file.path(path.out, paste(file.in,"output_best_visited_models.RData",sep="_"))
save(BestModels,file=ModelsWriteName)
return(BestModels)

}
###################################################################################
