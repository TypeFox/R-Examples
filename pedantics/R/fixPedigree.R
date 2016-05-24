`fixPedigree` <-
function(Ped, dat=NULL){

  if(is.null(dat)==FALSE&&is.null(dim(dat))==FALSE&&length(Ped[,1])!=length(dat[,1])) {
    cat(paste("Pedigree and cohorts differ in length.",'\n')); flush.console(); stop();
  }
  if(is.null(dat)==FALSE&&is.null(dim(dat))&&length(Ped[,1])!=length(dat)) {
    cat(paste("Pedigree and cohorts differ in length.",'\n')); flush.console(); stop();
  }

  names(Ped)<-c("id","dam","sire")
  ntotal<-length(Ped$id)*3
  IDs<-array(dim=ntotal)
  for(x in 1:length(Ped$id)) {
    IDs[x]<-as.character(Ped$id[x])
    IDs[x+ntotal]<-as.character(Ped$dam[x])
    IDs[x+ntotal*2]<-as.character(Ped$sire[x])
  }
  IDs<-as.data.frame(IDs)
  IDs<-unique(IDs)
  IDs<-subset(IDs,is.na(IDs)==FALSE)
  names(IDs)<-"id"
  IDs$dam<-Ped$dam[match(IDs$id,Ped$id)]
  IDs$sire<-Ped$sire[match(IDs$id,Ped$id)]
  orderPed<-function(ped){
    reorder<-ped[order(kindepth(ped[,1],ped[,2],ped[,3]), decreasing=FALSE),]
    return(reorder)
  }
  fixedPedigree<-orderPed(IDs)
  if(is.null(dat)==FALSE){
    if(names(dat)[1]=='id'|names(dat)[1]=='ID'|names(dat)[1]=='ids'|names(dat)[1]=='IDS'){
      for(x in 2:length(dat[1,])){
        fixedPedigree[,(3+x-1)]<-dat[match(fixedPedigree[,1],dat[,1]),x]
      }
    } else {
      cat(paste("No id column detected in dat, assuming same order as Ped.",'\n')); flush.console();
      dat$id<-Ped[,1]
      for(x in 1:(length(dat[1,])-1)){
        fixedPedigree[,(3+x-1)]<-dat[match(fixedPedigree[,1],dat$id),x]
      }
    }
  }
  for(x in 1:3) fixedPedigree[,x]<-as.factor(fixedPedigree[,x])
  fixedPedigree
}

