partigenD<-function(tr,grid=TRUE,zerorecs=FALSE)
{
d<-length(tr$N)
step<-stepcalc(tr$support,tr$N)

nodenum<-length(tr$left)
left<-tr$left
right<-tr$right
{
if (grid){
  low<-matrix(0,nodenum,d)
  upp<-matrix(0,nodenum,d)
  for (i in 1:nodenum){
    for (j in 1:d){
      low[i,j]<-tr$low[i,j]*step[j]+tr$support[2*j-1]
      upp[i,j]<-tr$upp[i,j]*step[j]+tr$support[2*j-1]
    }
  }
}
else{
  low<-tr$low
  upp<-tr$upp
}
}
mean<-tr$mean

ll<-leaflocs(left,right)
leafloc<-ll$leafloc
leafnum<-ll$leafnum

values<-matrix(0,leafnum,1)
recs<-matrix(0,leafnum,2*d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if ((mean[node]>0) || (zerorecs)){

     efek<-efek+1

     values[efek]<-mean[node]
 
     for (k in 1:d){

       recs[efek,2*k-1]<-low[node,k]
       recs[efek,2*k]<-upp[node,k]
     }

   }

   i<-i+1
}
values<-values[1:efek]
recs<-recs[1:efek,]
return(list(values=values,recs=recs))
}













