## this function is reading in the standard input of test set, and perform the optimization ##
optimize.features<-function(set,ranking,num=200,step=45,replication=5,list=FALSE){
  l1<-length(set[[2]])
  l2<-length(ranking)
  if(l1!=l2)
    return("alignment selection incorrect")
  data<-set[[2]][1:2,]
  Al<-matrix(ncol=length(data[1,]),nrow=2)
  for(i in 1:length(Al[1,])){
    Al[,i]<-as.character(set[[2]][c(3*ranking[i]+1,3*ranking[i]+2),i])
  } 

    M<-list()
    for(j in 1:replication){
      M[[j]]<-optimization.GA(Al=Al,data=data,num=num,step=step,plot=FALSE)
    }

  if(list==FALSE){
    para<-selection(M)
    return(para)
  }
  if(list==TRUE)
    return(M)
}