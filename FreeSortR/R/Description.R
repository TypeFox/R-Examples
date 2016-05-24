################################################################
#__________       Function  Description     _____________
#             Text description of a partition
################################################################
Description<-function(Partition,Labels=NULL){
  
  if (is.null(Labels)){
    Labels<-names(Partition)
  }
  nstimuli=length(Partition) 
  uniqPartition<-unique(Partition)
  nbr_groups=length(uniqPartition)
  
  # Generates the text representing the partition
  res=""
  for (group in (1:nbr_groups)){
    res=paste(res, "{",sep="")
    num<-which(Partition==uniqPartition[group])
    res=paste(res,Labels[num[1]],sep="")
    if (length(num)>1){
      for (i in (2:(length(num)))){
        res=paste(res, ", ",Labels[num[i]],sep="")
      }
    }  
    res=paste(res,"}",sep="")
  }  
  return(res)
}
