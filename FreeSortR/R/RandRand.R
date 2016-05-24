
###########################################################
# Computation of Rand Index
# between two partitions
###########################################################

RandRand<-function(Partition1,Partition2){
  
  
  nitem1<-length(Partition1)
  
  if (length(Partition2) != nitem1){
    print("problem of dimension of the two partitions")
  } else {
    #Partitions on common stimuli (case of incomplete sortings)
    common<-Partition1!=0 & Partition2!=0
    Partition1<-Partition1[common]
    Partition2<-Partition2[common]
    nitem<-sum(common)
    
    #Rand Classic
    t<-nitem*(nitem-1)/2

    A<-0  # agreements
    D<-0  # disagreements
    for (i in 1:(nitem-1)){
      for (j in (i+1):nitem){
        if ((Partition1[i]==Partition1[j]) & (Partition2[i]==Partition2[j])){
          A<-A+1
        } else {
          if ((Partition1[i]!=Partition1[j]) & (Partition2[i]!=Partition2[j])){
            D<-D+1
          }
        }
      }
    }
    Rand<-(A+D)/t
    return(Rand)
  }  
}
