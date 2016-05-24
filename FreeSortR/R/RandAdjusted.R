
###########################################################
# Computation of Adjusted Rand
# between two partitions
###########################################################

RandAdjusted<-function(Partition1,Partition2){
  
  
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
    Table<-table(Partition1,Partition2)
    nt<-rowSums(Table)
    pt<-colSums(Table)
    t<-nitem*(nitem-1)/2
    a<-sum( nt*(nt-1)/2) 
    b<-sum( pt*(pt-1)/2) 
    n<-sum( Table*(Table-1)/2)
    expectedrandindex<-a*b/t
    RandAdjusted<- (n-expectedrandindex)/(((a+b)/2)- expectedrandindex) 
    

    return(RandAdjusted)
  }  
}
