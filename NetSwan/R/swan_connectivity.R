swan_connectivity <-
function(g){
n<-length(V(g))
noeud<-V(g)
arc<-E(g)
fin <- matrix(ncol=1,nrow=n, 0)
dist<-distances(g)
dist[dist!=Inf]<-0
dist[dist==Inf]<-1
con<-sum(dist)
for(i in 1:n){
g2<-g
g2<-delete_vertices(g2, i)
dist2<-distances(g2)
dist2[dist2!=Inf]<-0
dist2[dist2==Inf]<-1
con2<-sum(dist2)
fin[i]<-con2-con
}
return(fin)
}
