swan_efficiency <-
function(g){
n<-length(V(g))
noeud<-V(g)
arc<-E(g)
fin <- matrix(ncol=1,nrow=n, 0)
dt=distances(g)
tot=sum(dt)
for(i in 1:n){
g2<-g
g2<-delete_vertices(g2, i)
dt2<-distances(g2)
tot2<-sum(dt2)
fin[i]<-tot2-(tot-sum(dt[i,])-sum(dt[,i]))
}
return(fin)
}
