partigen<-function(tr,grid=TRUE,zerorecs=FALSE)
{

d<-2
nodenum<-length(tr$left)
left<-tr$left
right<-tr$right
if (is.null(tr$step)) tr$step<-stepcalc(tr$support,tr$N)

if (grid){
  low<-matrix(0,nodenum,d)
  upp<-matrix(0,nodenum,d)
  for (i in 1:nodenum){
    for (j in 1:d){
      low[i,j]<-tr$low[i,j]*tr$step[j]+tr$support[2*j-1]
      upp[i,j]<-tr$upp[i,j]*tr$step[j]+tr$support[2*j-1]
    }
  }
}
else{
  low<-tr$low
  upp<-tr$upp
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
 
     recs[efek,1]<-low[node,1]
     recs[efek,2]<-upp[node,1]
  
     recs[efek,3]<-low[node,2]
     recs[efek,4]<-upp[node,2]
   }

   i<-i+1
}
values<-values[1:efek]
recs<-recs[1:efek,]
return(list(values=values,recs=recs))
}













