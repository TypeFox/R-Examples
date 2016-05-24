simutree<-function(tr,mcn,seedi)
{

set.seed(seedi)

step<-stepcalc(tr$support,tr$N)
d<-length(tr$N)
mcdendat<-matrix(0,mcn,d)

ll<-leaflocs(tr$left,tr$right)
leafnum<-ll$leafnum
leafloc<-ll$leafloc

p<-matrix(0,leafnum,1)
for (i in 1:leafnum){
  loc<-leafloc[i]
  p[i]<-tr$volume[loc]*tr$mean[loc]
}
p<-p/sum(p)

for (i in 1:mcn){
   ehto<-runif(1)
   alku<-0
   loppu<-p[1]
   lippu<-0
   for (j in 1:(leafnum-1)){
      if ((alku<=ehto) && (ehto<loppu)){
         loc<-leafloc[j]
         ran<-runif(d)
         for (k in 1:d){
            ala<-tr$suppo[2*k-1]+step[k]*tr$low[loc,k]
            yla<-tr$suppo[2*k-1]+step[k]*tr$upp[loc,k]
            ran[k]<-ala+(yla-ala)*ran[k]   
         }
         mcdendat[i,]<-ran
         lippu<-1
      }
      alku<-alku+p[j]
      loppu<-loppu+p[j+1]
   }      
   if (lippu==0){
         loc<-leafloc[leafnum]
         ran<-runif(d)
         for (k in 1:d){
            ala<-tr$suppo[2*k-1]+step[k]*tr$low[loc,k]
            yla<-tr$suppo[2*k-1]+step[k]*tr$upp[loc,k]
            ran[k]<-ala+(yla-ala)*ran[k]   
         }
         mcdendat[i,]<-ran
   }
}

return(mcdendat)

}

