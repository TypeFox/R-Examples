cvolumdya<-function(volofatom,component,AtomlistNext){
#
componum<-length(component)
volume<-matrix(0,componum,1)
#
# it is enough to calculate the number of aoms in each component
for (i in 1:componum){
   numofatoms<-0
   pointer<-component[i]
   while (pointer>0){
        numofatoms<-numofatoms+1
        pointer<-AtomlistNext[pointer]
   }
   volume[i]<-numofatoms*volofatom
}
return(volume)
}
