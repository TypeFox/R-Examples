dombin<-function(domOut){

nvar<-length(domOut$GD)

dombinout<-matrix(ncol=3,nrow=nvar*(nvar-1)/2)
colnames(dombinout)<-c("Comp","Cond","Gen")
dombinout<-data.frame(dombinout)
pair<-combn(1:nvar,2)
compDA<-domOut$DA[,c(-1,-2)]
compDA<-rbind(compDA,domOut$CD[1,])

for (i in 1:ncol(pair)){
  rownames(dombinout)[i]<-paste(names(domOut$GD)[pair[1,i]],names(domOut$GD)[pair[2,i]],sep=">")

  if (domOut$GD[pair[1,i]]>domOut$GD[pair[2,i]])
    dombinout$Gen[i]<-1
  else
    dombinout$Gen[i]<-0


  cdchk<-(domOut$CD[,pair[1,i]]>domOut$CD[,pair[2,i]])
  if (sum(cdchk) == nvar)
    dombinout$Cond[i]<-1
  else if (sum(cdchk) == 0)
    dombinout$Cond[i]<-0
  else
    dombinout$Cond[i]<-.5

  compDAl<-compDA[,c(pair[1,i],pair[2,i])]
  compDAl<-na.omit(compDAl)
  cdchk<-compDAl[,1]>compDAl[,2]
  if (sum(cdchk) == length(cdchk))
    dombinout$Comp[i]<-1
  else if (sum(cdchk) == 0)
    dombinout$Comp[i]<-0
  else
    dombinout$Comp[i]<-.5
}
return(dombinout)
}