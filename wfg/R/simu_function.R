network.simu <- function(nv=c(32,32,32,32), p.in=c(0.323,0.323,0.323,0.323), p.out=0.0625, p.del=0)
{
  n <- sum(nv) # number of vertices
  n.edge <- n*(n-1)/2
  A <- matrix(0, nrow=n, ncol=n) # initialize adjacency matrix
  l <- length(nv) # number of communities
  group <- c()

  for (k in 1:l) # obtain vector of membership
  {
    group <- c(group, rep(k,nv[k]))
  }

  from.mat <- replicate(n,group)
  to.mat <- t(from.mat)

  fromv <- from.mat[upper.tri(from.mat)]
  tov <- to.mat[upper.tri(to.mat)]

  index <- list()
  adj.vector <- rep(0, n.edge)

  for (i in 1:l)
  {
    index[[i]] <- which(tov==i & fromv==i)
    n.come <- length(index[[i]])
    unif.vector <- runif(n.come)
    com.adj.vector <- rep(0,n.come)
    com.adj.vector[unif.vector<p.in[i]] <- 1
    adj.vector[index[[i]]] <- com.adj.vector
  }

  index.btw <- which(tov!=fromv)

  n.btw <- length(index.btw)
  unif.vector <- runif(n.btw)
  btw.adj.vector <- rep(0,n.btw)
  btw.adj.vector[unif.vector<p.out] <- 1
  adj.vector[index.btw] <- btw.adj.vector

  # delete missing links
  ne <- sum(adj.vector) # number of edges
  nd <- round(ne*p.del) # number of edges to be deleted
  index.del <- sample(1:ne, nd) # randomly pick nd links to delete
  index.edge <- which(adj.vector==1)
  adj.vector[index.edge[index.del]] <- 0

  A[upper.tri(A)] <- adj.vector
  A <- A+t(A)

  net <- graph.adjacency(A, mode="undirected")

  value <- list()
  value$net <- net
  value$group <- group
  return(value)
}
