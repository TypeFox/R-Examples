drawcart<-function(treeseq,lehtilkm,suppo,plkm){
#Makes data for drawing a perspective plot.
#
#plkm on kuvaajan hilan pisteiden lkm
#
#koe<-drawcart(dendat,lehtilkm=3,treeseq,epsi=0,plkm=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60) 
#
# Haetaan ensin haluttu puu puitten jonosta
tree<-treeseq$tree
delnodes<-treeseq$delnodes
delend<-treeseq$delend
leafs<-treeseq$leafs
indeksi<-detsi(leafs,lehtilkm)
endi<-delend[indeksi]
if (endi>0){  #if there is something to remove
  indeksit<-delnodes[1:endi]
  re<-remnodes(tree$left,tree$right,indeksit)  
  tree$left<-re$left
  tree$right<-re$right
}
# Tehdaan binaaripuusta paloittain vakio
pv<-partition(tree,suppo)       
recs<-pv$recs
values<-pv$values
#
ans<-drawgene(values,recs,plkm)
return(list(x=ans$x,y=ans$y,z=ans$z))
}










