random.attack <- function(adj.mat,max.remove,nsim)
{
  n.nodes<-dim(adj.mat)[1]

  size.sim<-rep(0,n.nodes)
  mean.size.large.connex<-rep(0,n.nodes)
  mean.charac.path.length<-rep(0,n.nodes)

  for(i in 1:nsim){

    remove<-rep(0,n.nodes)
    rem.nodes<-rep(0,n.nodes)
    charac.path.length<-rep(0,n.nodes)
    size.large.connex<-rep(0,n.nodes)

    for(j in 1:max.remove)
    {
      k<-0
      tmp<-ceiling(runif(1,0,(n.nodes-j+1)))
      while(tmp>0){
        k<-k+1
        if(remove[k]==0) tmp<-tmp-1
      }
      remove[k]<-1
      rem.nodes[j]<-k
      rm<-node.attack(adj.mat,k)
      adj.mat<-rm$new.mat

      size.large.connex[j]<-rm$size.large.connex
      charac.path.length[j]<-rm$charac.path.length

      if(charac.path.length[j]!=-1)
      {
        mean.charac.path.length[j]<-mean.charac.path.length[j]+charac.path.length[j]
        size.sim[j]<- size.sim[j]+1
      }
    }
    mean.size.large.connex<-mean.size.large.connex+size.large.connex
  }
  mean.charac.path.length<-mean.charac.path.length/ size.sim
  mean.size.large.connex<-mean.size.large.connex/nsim
  list(size.large.connex=size.large.connex,charac.path.length=charac.path.length,rem.nodes=rem.nodes,mean.size.large.connex=mean.size.large.connex,mean.charac.path.length=mean.charac.path.length)
}
















