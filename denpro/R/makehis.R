makehis<-function(regdat)
{
xlkm<-length(regdat$hila[,1])  #muuttujien lkm
valipit<-matrix(0,1,xlkm)
i<-1
while (i<=xlkm){
  if (regdat$hila[i,1]>1) 
    valipit[i]<-(regdat$hila[i,3]-regdat$hila[i,2])/(regdat$hila[i,1]-1)
  i<-i+1
}
lnum<-length(regdat$ind[,1])   #length(regdat$dep)
items<-matrix(0,lnum,2*xlkm)
arvot<-matrix(0,lnum,1)
i<-1
while (i<=lnum){
  arvot[i]<-regdat$dep[i]  
  j<-1
  while (j<=xlkm){
    items[i,2*j-1]<-regdat$ind[i,j]-valipit[j]/2
    items[i,2*j]<-regdat$ind[i,j]+valipit[j]/2
    j<-j+1
  }
  i<-i+1
}
return(list(values=arvot,recs=items))
}


