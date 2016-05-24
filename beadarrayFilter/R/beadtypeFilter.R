beadtypeFilter <-
function(beadsum,Quantile=1,keepData=TRUE,delta=0.5)
{
allo<-c("ExpressionSetIllumina","LumiBatch", "data.frame")
checkclass <- class(beadsum)
if(!(checkclass %in% allo)) {stop("beadsum Object must be an ExpressionSetIllumina or a data.frame ")}
else 
if (checkclass =="ExpressionSetIllumina"){  
###delete rows with missing data.
aaa<-na.omit(data.frame(I(rownames(exprs(beadsum))),exprs(beadsum)))
    eSet <- na.omit(exprs(beadsum))
    ##obtain the std dev and convert them to stderr
    stdev <- na.omit(se.exprs(beadsum))
    nSet <- na.omit(attributes(beadsum)$assayData$nObservations)
    ###calculate stderr from the std Dev stderr=std Dev/sqrt(nset)
    seSet<-stdev/sqrt(nSet)
    ProbeID <- aaa[,1]
    #group <- c(1:dim(eSet)[2])
iccResults<-iccFun(eSet=eSet,seSet,nSet=nSet,ProbeID =ProbeID,iccQuant=Quantile,diffIcc=FALSE,keepData=FALSE)
 informID<- subset(iccResults$icc, iccResults$icc[,2]>=delta)[,1]
x<-beadsum[informID, ]
}
else if (checkclass =="LumiBatch"){ 

eSet <- exprs(beadsum)
seSet<- se.exprs(beadsum) ##this data has stderr and not stdev
nSet <- beadNum(beadsum)
 #group <- c(1:dim(eSet)[2])
ProbeID =fData(beadsum)$ProbeID 
iccResults<-iccFun(eSet=eSet,seSet,nSet=nSet,ProbeID =ProbeID,iccQuant=Quantile,diffIcc=FALSE,keepData=FALSE)
 informID<- subset(iccResults$icc, iccResults$icc[,2]>=delta)[,1]
x<-beadsum[informID, ]
}
else if (checkclass =="data.frame"){  
####low level function that uses a matrix instead of an expressionSet
eSet <- beadsum[, grep("Signal", names(beadsum))]
    seSet <- beadsum[, grep("STDERR", names(beadsum))]
    nSet <- beadsum[, grep("NBEADS", names(beadsum))]
    ProbeID <- beadsum[, 1]
iccResults<-iccFun(eSet=eSet,seSet=seSet,nSet=nSet,ProbeID =ProbeID,iccQuant=Quantile,diffIcc=FALSE,keepData=FALSE)
 informID<- subset(iccResults$icc, iccResults$icc[,2]>=delta)[,1]
x<-subset(beadsum,beadsum[,1] %in% informID)
}
if (keepData)
    return(list(InformProbeNames=informID,informData=x))
else  return(list(InformProbeNames=informID,informData=NULL))
}
