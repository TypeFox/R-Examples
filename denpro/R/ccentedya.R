ccentedya<-function(volofatom,component,AtomlistNext,AtomlistAtom,
volume,minim,h,delta,index,d){
#
componum<-length(component)
center<-matrix(0,componum,d)
#
for (i in 1:componum){
   curcente<-0
   pointer<-component[i]
   while (pointer>0){
        atompointer<-AtomlistAtom[pointer]
        inde<-index[atompointer,]
        newcente<-minim-h+delta*inde
        curcente<-curcente+newcente
        pointer<-AtomlistNext[pointer]
   }
   center[i,]<-volofatom*curcente/volume[i]
}
return(t(center))
}
