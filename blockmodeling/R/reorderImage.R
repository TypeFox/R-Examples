reorderImage<-function(IM,oldClu,newClu){
if(crand2(oldClu,newClu)!=1)stop("Old and new clu's are not compatibale (crand index is not 1)!\n")
newOrder<-which(table(oldClu,newClu)>0,arr.ind=TRUE)[,1]
return(IM[newOrder,newOrder])
}
