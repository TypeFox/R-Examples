varScore <- function(Data, Forward, Reverse=NULL, Range = NULL, Prorate = TRUE, MaxMiss = .20)
{
  #select relevant items
  d = Data[,c(Forward, Reverse)]
  
  #check for out of range
  if (!is.null(Range)){
    if (min(d, na.rm=TRUE) < Range[1] || max(d, na.rm=TRUE) > Range[2]){
      stop('Item score(s) out of range')
    }
  }
  
  #check that length of Range == 2 if Reverse is not null
  if (!is.null(Reverse) && length(Range) !=2) {
    stop('Must specify item range (Range) to reverse score items')
  }
  
  #Reverse score relevant items
  if (!is.null(Reverse)){
    for (v in Reverse) {
      d[,v] = (Range[1] + Range[2]) - d[,v]
    }   
  }
  
  if (Prorate){
    Total = rowMeans(d, na.rm=TRUE)*dim(d)[2]
  }
  else{
    Total = rowSums(d, na.rm=TRUE)
  }
  
  #count missing and set > MaxMiss to NA
  MissCount = rowSums(is.na(d))
  MissCount = MissCount/dim(d)[2]
  Total[MissCount > MaxMiss] = NA
  
  return(Total)
}