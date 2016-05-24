
###########################################################
# Computation of Rand and Adjusted Rand
# between two partitions
###########################################################

RandIndex<-function(Partition1,Partition2){
  
  
  nitem1<-length(Partition1)
  
  if (length(Partition2) != nitem1){
    print("problem of dimension of the two partitions")
  } else {
    #Partitions on common stimuli (case of incomplete sortings)
    common<-Partition1!=0 & Partition2!=0
    Partition1<-Partition1[common]
    Partition2<-Partition2[common]
    nitem<-sum(common)
    
    #Rand Adjusted

    RandAdjusted<-RandAdjusted(Partition1,Partition2)
    
    # Rand
    
    Rand<-RandRand(Partition1,Partition2)
    
    return(list(Rand=Rand,AdjustedRand=RandAdjusted))
  }  
}
