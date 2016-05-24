ABCRemoveSmallYields=function(Data,CumSumSmallestPercentage=0.5){
# res = ABCRemoveSmallYields(Data,CumSumSmallestPercentage)
# Data cleanning for ABC analysis
# the smallest data up to a cumulated sum of less than CumSumSmallestPercentage percent
# the total sum (yield) is removed 
# negative data values and NaN are treated as zeros
#
# INPUT
# Data(1:n)           the data set, may contain NaN, negative values and very small values
# 
# OPTIONAL
# CumSumSmallestPercentage  (default =0.5),the smallest data up to a cumulated sum of less than CumSumSmallestPercentage
#
# OUTPUT
# SubstantialData(1:n1)         columnvector containing Data>=0 and zeros for all NaN and negative values in Data(1:n)  
# Data2SubstantialInd           Index such that SubstantialData = nantozero(Data(Data2SubstantialInd));
# RemovedInd                    Data(RemovedInd) is the data that has been removed 
# author: MT 08/2015
   CleanDatares=ABCcleanData(Data)
   CleanData=CleanDatares$CleanedData
  #die kleinsten Daten, die zusammen weniger als CumSumSmallestPercentage# ausmachen identifizieren
  SortedData=sort(CleanData,decreasing=FALSE)
  TotalYield = sum(SortedData)
  CumSumPercentage=round(cumsum(SortedData/TotalYield*100),0)
 SmallInd=which(CumSumPercentage<CumSumSmallestPercentage) #In Prozent
# Wenn Welche gefunden werden
if(length( SmallInd) >0){
     # print('Removing the smallest data up to a cumulated sum of less than 0.5% of the total sum (yield):')
      SchwellenIndex=tail(SmallInd,1)+1 # Index des letzten zu kleinen Datensatzes
      Schwelle = SortedData[SchwellenIndex]      # der Datensatz der gerade noch verbleiben kann
      Data2CleanInd = which(CleanData>=Schwelle)  # diese Daten koennen bleiben
      RemovedInd    = which(CleanData<Schwelle)  # diese Daten fallen weg
      CleanData   = CleanData[Data2CleanInd]  # bereinigen      
     # print(paste0(length(RemovedInd),' items removed.'))
}else{
      Data2CleanInd=CleanDatares$Data2CleanInd
      RemovedInd=setdiff(CleanDatares$Data2CleanInd,1:length(Data))
}
 
return(list(SubstantialData=CleanData,Data2SubstantialInd=Data2CleanInd,RemovedInd=RemovedInd))
}
  