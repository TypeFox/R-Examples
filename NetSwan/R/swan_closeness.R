swan_closeness <-
function(g){
n<-length(V(g))
noeud<-V(g)
arc<-E(g)
fin <- matrix(ncol=1,nrow=n, 0)
clos<-distances(g)
closb<-1/clos
closb[closb==Inf]<-0
tot<-sum(closb)
for(i in 1:n){
g2<-g
g2<-delete_vertices(g2, i)
clos2<-distances(g2)
closb2<-1/clos2
closb2[closb2==Inf]<-0
tot2<-sum(closb2)
fin[i]<-tot2-(tot-sum(closb[i,])-sum(closb[,i]))
}
return(fin)
}
