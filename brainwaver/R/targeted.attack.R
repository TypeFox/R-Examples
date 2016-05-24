targeted.attack <- function(adj.mat,max.remove)
{
  n.nodes<-dim(adj.mat)[1]

  size.sim<-rep(0,n.nodes)
  remove<-rep(0,n.nodes)
  rem.nodes<-rep(0,n.nodes)

  charac.path.length<-rep(0,n.nodes)
  size.large.connex<-rep(0,n.nodes)

  in.degree<-rowSums(adj.mat)
  in.degree.sort<-sort(in.degree,index.return=TRUE,decreasing=TRUE)

  for(j in 1:max.remove)
  {

    tmp<-in.degree.sort$ix[j]

    rem.nodes[j]<-tmp
    k<-tmp

    rm<-node.attack(adj.mat,k)
    adj.mat<-rm$new.mat
    size.large.connex[j]<-rm$size.large.connex
    charac.path.length[j]<-rm$charac.path.length

  }

  list(size.large.connex=size.large.connex,charac.path.length=charac.path.length,rem.nodes=rem.nodes)
}





















