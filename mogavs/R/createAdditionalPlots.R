createAdditionalPlots <-
function(mogavs, epsilonBand=0, kBest=1, method=c("MSE","kBest")){
  method<-tolower(method)
  archiveSet<-mogavs$archiveSet
  obj1ArchiveSet<-mogavs$obj1ArchiveSet
  obj2ArchiveSet<-mogavs$obj2ArchiveSet
  N<-ncol(mogavs$archiveSet)
  
  if(method=="mse" && missing(epsilonBand)){
    warning("Arg epsilonBand not supplied, defaulting to zero",call.=FALSE)
    epsilonBand<-0
  }
  if(method=="kbest" && missing(kBest)){
    warning("Arg kBest not supplied, defaulting to one")
    kBest<-1
  }
  sizeArchive<-nrow(archiveSet)
  epsMembers<-list()
  epsMembersTemp<-list()
  minMSE<-list(NA)
  
  for(i in 1:N){
    #Finds all members with i number of variables
    iMembers<-which(obj1ArchiveSet==i)
    if(length(iMembers)==0){
      break
    }
    #Finds the member with i number of variables and minimum MSE
    minMSE[[i]]<-min(obj2ArchiveSet[iMembers])
    #Finds the tag of all members with i variables within epsilon range
    if(length((minMSE[[i]]))>0){
      epsMembersTemp[[i]]<-obj2ArchiveSet[iMembers]-minMSE[[i]]<=epsilonBand
      #Stores all the members with i number of variables within epsilon range
      epsMembers[[i]]<-iMembers[epsMembersTemp[[i]]]
      } else{
        epsMembersTemp[[i]]<-""
        epsMembers[[i]]<-""
      }
  }
  if(method=="mse"){
  temp<-data.frame()
  for(i in 1:N){
    if(length(epsMembers[[i]])>0){
      temp<-rbind(temp,cbind(obj1ArchiveSet[epsMembers[[i]]],obj2ArchiveSet[epsMembers[[i]]]))
    }
    
  }
  plot(obj1ArchiveSet,obj2ArchiveSet,col="red",pch=8)
  points(temp,col="blue",pch=8)
  legend("topright",c("Entire Archive Set",paste("Members in Epsilon (=",epsilonBand,") band",sep="")),col=c("red","blue"),pch=c(8,8))
  rm(temp)
}
if(method=="kbest"){
  #Create new figure to plot k best models for each number of variables
  plot(obj1ArchiveSet,obj2ArchiveSet,col="red",pch=8)
  orderedMSE<-list()
  kBestMembersTemp<-list()
  kBestMembers<-list()
  for(i in 1:N){
    #Finds all members with i number of variables
    iMembers<-which(obj1ArchiveSet==i)
    #Finds the member with i number of variables and ordered MSE
    orderedMSE[[i]]<-sort(obj2ArchiveSet[iMembers])
    I<-order(obj2ArchiveSet[iMembers])
    #finds the tag of all members that are kBest
    if(length(orderedMSE[[i]])>0){
      kBestMembersTemp[[i]]<-iMembers[I]
      s<-min(length(kBestMembersTemp[[i]]),kBest)
      kBestMembers[[i]]<-kBestMembersTemp[[i]][1:s]
    }
    else {
      kBestMembers[[i]]<-"NA"
    }
  }
  
  for(i in 1:N){
    if(length(kBestMembers[[i]])>0){
      #Plot the members with i number of variables and kBest
      points(obj1ArchiveSet[kBestMembers[[i]]],obj2ArchiveSet[kBestMembers[[i]]],col="blue",pch=8)
      legend("topright",c("Entire Archive Set",paste("kBest(=",kBest,") members for each level of variables",sep="")),col=c("red","blue"),pch=c(8,8))
      
    }
  }
}

}
