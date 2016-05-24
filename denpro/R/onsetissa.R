onsetissa<-function(kandi,h,delta,minim,
brnode,component,
index,
AtomlistAtom,AtomlistNext){
#
itis<-F
d<-length(minim)
# 
node<-brnode
compo<-component[node]
ato<-compo                          #ato is pointer to "value"
while ((ato>0) && !(itis)){
    inde<-index[AtomlistAtom[ato]]
    keski<-minim-h+delta*inde
    for (din in 1:d){
      if ((kandi[din]>=(keski[din]-delta[din]/2)) &&   
          (kandi[din]<=(keski[din]+delta[din]/2))){
               itis<-T
      }
    }
    ato<-AtomlistNext[ato]
}
#
return(itis)
}
