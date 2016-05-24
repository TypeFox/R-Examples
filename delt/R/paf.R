paf<-function(a,b){
#
rown<-dim(a)[1]
coln<-dim(a)[2]
res<-matrix(0,rown,coln)
#
counter<-0
for (i in 1:rown){
  if (b[i]==1){
     counter<-counter+1 
     res[counter,]<-a[i,]
  }
}
#
res<-res[1:counter,]
return(res)
}
