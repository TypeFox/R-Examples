ABCcleanData=function(Data){
# V= ABCcleanData(Data)
# Data cleanning for ABC analysis
# only the first column of Data is used
# Data <0 are set to zero,  NA in Data are set to zero
# if RemoveSmallYields ==TRUE => the smallest data up to a cumulated sum of less than
# 0.5# of the total sum (yield) is removed 
#
# INPUT
# Data(1:n)           the data set, may contain NaN, negative values and very small values
# 
# OUTPUT      List V with
# V$CleanedData(1:n1)   columnvector containing Data>=0 and zeros for all NaN and negative values in Data(1:n)  
# V$Data2CleanInd       Index such that CleanedData = nantozero(Data(Data2CleanInd));
#
# author MT 01/2015, reimplemented from ALUs matlab version
# 1.Editor: MT 08/2015 RemoveSmallYields in eigene Funktion verschoben
if(!is.vector(Data)){
  n=nrow(Data)
  d=ncol(Data)
  warning('Only vectors should be used!')
  if(d>1){ #Data is Matrix or data.frame
    warning('Using only first column of data')
    UncleanData=as.vector(Data[,1]) # use only first column
  }else{
    UncleanData=Data
  }
}else{
  UncleanData=Data
}
  UncleanData=as.numeric(unname(UncleanData))#Automatische Umwandlung voh Chars/string in NA
  rowsbefore=length(UncleanData)
#NA durch Nullen ersetzen
  nabools=is.finite(UncleanData)
  Data2CleanInd=which(nabools==FALSE)
  CleanData=UncleanData
  if(length(Data2CleanInd)) CleanData[Data2CleanInd]=0
#Negative Werte entfernen

  DataNeg=CleanData[CleanData<0]
  cols=1
  bools=CleanData %in% DataNeg
  CleanData[bools]<-0

  rows=rowsbefore-sum(bools)-sum(!nabools)
if(rowsbefore>rows){
     warning(paste0(rows,' of ',rowsbefore,' items are positive and beeing used for further calculations.'))
    # warning('Please use Data[Data>0], before using Data[Aind] etc.')    
}

return(list(CleanedData=CleanData,Data2CleanInd=Data2CleanInd))
}
