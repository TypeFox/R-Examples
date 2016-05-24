#given an adjacency matrix, remove the specified node,
#and then return the new matrix, together with size of its largest connected component
#and the new characteristic path length
node.attack <- function(adj.mat,node.sup)
{
  n.nodes<-dim(adj.mat)[1]

  adj.mat[node.sup,]<-rep(0,n.nodes)
  adj.mat[,node.sup]<-rep(0,n.nodes)

  z<-.C("Rconnex_components",
        as.integer(n.nodes),
        as.integer(adj.mat),
        acc=double(n.nodes+3), PACKAGE="brainwaver")

  size.large.connex<-z$acc[n.nodes+2]
  idx.large.connex<-z$acc[n.nodes+1]

  z<-.C("Rcharac_path_length",
        as.integer(idx.large.connex),   
        as.integer(n.nodes),
        as.integer(size.large.connex),
        as.integer(adj.mat),
        Lp=double(1), PACKAGE="brainwaver")

  charac.path.length<-z$Lp
  list(new.mat=adj.mat, size.large.connex=size.large.connex,charac.path.length=charac.path.length)
}



















