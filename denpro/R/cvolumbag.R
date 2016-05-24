cvolumbag<-function(component,AtomlistAtom,AtomlistNext,low,upp,steppi)
{
d<-dim(low)[2]

componum<-length(component)
volume<-matrix(0,componum,1)

for (i in 1:componum){
   curvolu<-0
   pointer<-component[i]
   while (pointer>0){
        atto<-AtomlistAtom[pointer]

        vol<-1
        for (j in 1:d){
            vol<-vol*(upp[atto,j]-low[atto,j])*steppi[j]
        }

        curvolu<-curvolu+vol
        pointer<-AtomlistNext[pointer]
   }
   volume[i]<-curvolu
}
return(volume)
}
