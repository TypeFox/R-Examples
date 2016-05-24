# ' compare 0 values in two vectors
# ' @param trueA first vector
# ' @param Aalgo second vector
# ' @param taux boolean. Computes the ratio of each statistic or not.
compare_zero<-function(trueA=trueA,Aalgo=Aalgo,taux=FALSE){
  quivrai0=which(trueA==0)
  if(length(quivrai0)>0){
    nbbon0=length(which(Aalgo[quivrai0]==0))
    nbbon1=length(which(Aalgo[-quivrai0]!=0))
    nbfaux0=length(which(Aalgo[-quivrai0]==0))
    nb0mank=length(quivrai0)-nbbon0
  }else{
    nbbon0=0
    nbbon1=length(which(Aalgo!=0))
    nbfaux0=length(which(Aalgo==0))
    nb0mank=0
  }
  if(taux){
    tauxbon0=nbbon0/length(quivrai0)
    tauxfaux0=nbfaux0/length(which(Aalgo==0))
    return(list(true0=nbbon0,true1=nbbon1,false0=nbfaux0,false1=nb0mank,ratio_true0=tauxbon0,ratio_false0=tauxfaux0))
  }else{
    return(list(true0=nbbon0,true1=nbbon1,false0=nbfaux0,false1=nb0mank))    
  }
}