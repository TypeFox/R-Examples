partition<-function(et,grid=TRUE,zerorecs=FALSE)
{

d<-length(et$support)/2
nodenum<-length(et$left)
left<-et$left
right<-et$right
if (is.null(et$step)) et$step<-stepcalc(et$support,et$N)

if (grid){
  low<-matrix(0,nodenum,d)
  upp<-matrix(0,nodenum,d)
  for (i in 1:nodenum){
    for (j in 1:d){
      low[i,j]<-et$low[i,j]*et$step[j]+et$support[2*j-1]
      upp[i,j]<-et$upp[i,j]*et$step[j]+et$support[2*j-1]
    }
  }
}
else{
  low<-et$low
  upp<-et$upp
}

mean<-et$mean

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
     
     for (dd in 1:d){
         recs[efek,2*dd-1]<-low[node,dd]
         recs[efek,2*dd]<-upp[node,dd]
     }
   }

   i<-i+1
}

values<-values[1:efek]
if (efek==1) recs<-matrix(recs[1:efek,],1,2*d)
else recs<-recs[1:efek,]

return(list(values=values,recs=recs,support=et$support))
}













