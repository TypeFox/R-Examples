ccentebag<-function(component,AtomlistAtom,AtomlistNext,low,upp,volume,
step,suppo)
{
d<-dim(low)[2]

componum<-length(component)
center<-matrix(0,componum,d)

for (i in 1:componum){
   curcente<-matrix(0,d,1)
   pointer<-component[i]
   while (pointer>0){
        atompointer<-AtomlistAtom[pointer]
        
        newcente<-matrix(0,d,1)
        for (j in 1:d){
            # calculate 1st volume of d-1 dimensional rectangle where
            # we have removed j:th dimension

            vol<-1
            k<-1
            while (k<=d){
               if (k!=j){
                  vol<-vol*(upp[atompointer,k]-low[atompointer,k])*step[k]
               }
               k<-k+1
            }

            ala<-suppo[2*j-1]+step[j]*low[atompointer,j]
            yla<-suppo[2*j-1]+step[j]*upp[atompointer,j]
            newcente[j]<-vol*(yla^2-ala^2)/2
        }

        curcente<-curcente+newcente
        pointer<-AtomlistNext[pointer]
   }
   center[i,]<-curcente/volume[i]
}
return(t(center))
}

