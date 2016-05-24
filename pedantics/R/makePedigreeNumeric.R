`makePedigreeNumeric` <-
function(id,sire,dam,missingVal=NULL){
  p<-cbind(id,sire,dam)
  pf<-as.factor(p)
  pn<-as.numeric(pf)
  for(i in 1:length(pn)) {if(is.null(missingVal)==FALSE&is.na(pn[i])==TRUE) pn[i]<- missingVal}
  k<-as.data.frame(cbind(as.numeric(pn),as.character(pf)))
  k<-unique(k);  k[,1]<-as.numeric(as.character(k[,1]))
  names(k)<-c("pn","pf")
  ped<-as.data.frame(matrix(pn,length(id),3,byrow=FALSE))
  names(ped)<-c("id","sire","dam")
  r<-list(numericPedigree=ped,idKey=k)
  r
}

