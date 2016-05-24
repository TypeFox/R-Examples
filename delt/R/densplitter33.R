densplitter33<-function(dendat,minobs=1,method="loglik")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

suppo<-matrix(0,2*d,1)  
indet<-matrix(0,2*d,1)  
for (i in 1:d){
    suppo[2*i-1]<-min(dendat[,i])   
    suppo[2*i]<-max(dendat[,i])
    indet[2*i-1]<-seq(1,n)[(suppo[2*i-1]==dendat[,i])]
    indet[2*i]<-seq(1,n)[(suppo[2*i]==dendat[,i])]
}
notindet<-setdiff(seq(1,n),indet)
inside<-dendat[notindet,]
m<-dim(inside)[1]

x<-matrix(0,m+1,d)  # x contains split points and boundaries
for (i in 1:d){
    ordi<-order(inside[,i])
    x[1,i]<-min(dendat[,i])   
    x[m+1,i]<-max(dendat[,i])
    ala<-inside[ordi[1:(m-1)],i]
    yla<-inside[ordi[2:m],i]
    x[2:m,i]<-(ala+yla)/2
}

maksi<-2*n  #n^2
currecs<-matrix(0,maksi,2*d)
for (i in 1:d){
   currecs[1,2*i-1]<-1
   currecs[1,2*i]<-m+1
}
pinin<-1
finrecs<-matrix(0,maksi,2*d)
saatu<-0

while (pinin>0){
   rec<-currecs[pinin,]
   pinin<-pinin-1
   
   fs<-findsplitter(x,rec,n,method)     
   direc<-fs$vec
   point<-fs$val
   leftrec<-rec
   leftrec[2*direc]<-point
   rightrec<-rec
   rightrec[2*direc-1]<-point 

   leftdiffes<-matrix(0,d,1)
   rightdiffes<-matrix(0,d,1)
   for (dd in 1:d){
       leftdiffes[dd]<-leftrec[2*dd]-leftrec[2*dd-1]
       rightdiffes[dd]<-rightrec[2*dd]-rightrec[2*dd-1]
   }
   lkmleft<-min(leftdiffes)
   lkmright<-min(rightdiffes)
   if ((lkmleft>minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-leftrec
      currecs[pinin+2,]<-rightrec
      pinin<-pinin+2
   }
   if ((lkmleft>minobs)&&(lkmright<=minobs)){
      currecs[pinin+1,]<-leftrec
      pinin<-pinin+1
      finrecs[saatu+1,]<-rightrec
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-rightrec
      pinin<-pinin+1
      finrecs[saatu+1,]<-leftrec
      saatu<-saatu+1
   }
}

recs<-finrecs[1:saatu,]
if (saatu==1) recs<-matrix(recs,1,2*d)
truerecs<-recs
for (dd in 1:d){ 
    truerecs[,2*dd-1]<-x[recs[,2*dd-1],dd]
    truerecs[,2*dd]<-x[recs[,2*dd],dd]
}

return(list(recs=truerecs,support=suppo))
}





