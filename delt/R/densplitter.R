densplitter<-function(dendat,minobs=1,method="loglik",dyadic=FALSE)
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

grid<-matrix(0,m+1,d)  # grid contains split points and boundaries
for (i in 1:d){
    ordi<-order(inside[,i])
    grid[1,i]<-min(dendat[,i])   
    grid[m+1,i]<-max(dendat[,i])
    ala<-inside[ordi[1:(m-1)],i]
    yla<-inside[ordi[2:m],i]
    grid[2:m,i]<-(ala+yla)/2
}

maksi<-2*n  #n^2
currecs<-matrix(0,maksi,2*d)
currecs[1,]<-suppo
pinin<-1
finrecs<-matrix(0,maksi,2*d)
saatu<-0
curobs<-matrix(FALSE,maksi,m)
curobs[1,]<-TRUE
curdown<-matrix(0,maksi,d)
curhigh<-matrix(0,maksi,d)
curdown[1,]<-rep(1,d)
curhigh[1,]<-rep(m+1,d)
findown<-matrix(0,maksi,d)
finhigh<-matrix(0,maksi,d)

while (pinin>0){
   rec<-currecs[pinin,]   
   obs<-curobs[pinin,]
   recdown<-curdown[pinin,]
   rechigh<-curhigh[pinin,]
   pinin<-pinin-1
     
   x<-inside[obs,]
   if (dyadic) fs<-findsplitter.dyadic(grid,x,rec,n,method,minobs,recdown,rechigh)     
   else fs<-findsplitter(grid,x,rec,n,method,minobs,recdown,rechigh)     
   direc<-fs$vec
   point<-fs$val
   gridpoint<-fs$valio

   leftobs<-(obs&(inside[,direc]<point))
   rightobs<-(obs&(inside[,direc]>point))

   leftrec<-rec
   leftrec[2*direc]<-point
   rightrec<-rec
   rightrec[2*direc-1]<-point 

   leftdown<-recdown
   lefthigh<-rechigh
   lefthigh[direc]<-gridpoint

   rightdown<-recdown
   righthigh<-rechigh
   rightdown[direc]<-gridpoint

   lkmleft<-sum(leftobs)
   lkmright<-sum(rightobs)
   if ((lkmleft>minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-leftrec
      currecs[pinin+2,]<-rightrec
      curobs[pinin+1,]<-leftobs
      curobs[pinin+2,]<-rightobs
      curdown[pinin+1,]<-leftdown
      curhigh[pinin+1,]<-lefthigh
      curdown[pinin+2,]<-rightdown
      curhigh[pinin+2,]<-righthigh
      pinin<-pinin+2
   }
   if ((lkmleft>minobs)&&(lkmright<=minobs)){
      currecs[pinin+1,]<-leftrec
      curobs[pinin+1,]<-leftobs
      curdown[pinin+1,]<-leftdown
      curhigh[pinin+1,]<-lefthigh
      pinin<-pinin+1
      finrecs[saatu+1,]<-rightrec
      findown[saatu+1,]<-rightdown
      finhigh[saatu+1,]<-righthigh
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-rightrec
      curobs[pinin+1,]<-rightobs
      curdown[pinin+1,]<-rightdown
      curhigh[pinin+1,]<-righthigh
      pinin<-pinin+1
      finrecs[saatu+1,]<-leftrec
      findown[saatu+1,]<-leftdown
      finhigh[saatu+1,]<-lefthigh
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright<=minobs)){
      finrecs[saatu+1,]<-leftrec
      finrecs[saatu+2,]<-rightrec
      findown[saatu+1,]<-leftdown
      findown[saatu+2,]<-rightdown
      finhigh[saatu+1,]<-lefthigh
      finhigh[saatu+2,]<-righthigh
      saatu<-saatu+2
   }
}

recs<-finrecs[1:saatu,]
down<-findown[1:saatu,]
high<-finhigh[1:saatu,]
if (saatu==1){
   recs<-matrix(recs,1,2*d)
   down<-matrix(down,1,d)
   high<-matrix(high,1,d)
}

return(list(recs=recs,support=suppo,grid=grid,down=down,high=high))
}





