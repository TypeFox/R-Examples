`makePedigreeFactor` <-
function(id,sire,dam,key){
  p<-as.data.frame(as.factor(cbind(id,sire,dam)))
  p$ids<-key$pf[match(p[,1],key$pn)]
  ped<-as.data.frame(matrix(p$ids,length(id),3,byrow=FALSE))
  for(x in 1:3) ped[,x]<-as.factor(ped[,x])
  names(ped)<-c("id","sire","dam")
  ped
}

