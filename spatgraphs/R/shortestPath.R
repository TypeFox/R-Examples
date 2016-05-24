#' shortest path on the graph
#'
#' Djikstra's algorithm
#'
#' @param i index from
#' @param j index to
#' @param g sg object
#' @param x optional point pattern from which g was computed
#' @param dbg verbose
#'
#' @export

shortestPath<-function(i, j, g, x=NULL, dbg=FALSE)
{
  if(missing(g))stop("Give a graph g.")
  if(!(i%in%1:g$N)| !(j%in%1:g$N) | i==j) stop("Give i,j different and between 1,...,n.")
  g<-sg2sym(g)
  e<-spatcluster(g)

  for(k in 1:length(e$clusters))
  {
    if(i%in%e$clusters[[k]]) break;
  }

  if(!(j%in%e$clusters[[k]])) return(list(d=Inf, path=NA));

  if(is.null(x)) d<-matrix(1, g$N, g$N) - diag(g$N)
  else d <- as.matrix(dist(sg_parse_coordinates(x),upper=TRUE))

  cluster<-e$clusters[[k]]

  first<-g$edges[[i]]

  if(j %in%first) return(list(d=d[i,j],path=c(i,j)))


  #loop<-TRUE

  dists<-rep(Inf, length(cluster))
  previous<-rep(NA,length(cluster))
  ii<-which(i==cluster)
  dists[ii]<-0
  Q<-cluster
  left<-rep(TRUE,length(Q))

  while(sum(left)>0)
  {
    u<-which(min(dists[left])==dists)
    u<-u[ceiling(runif(1)*length(u))]
    uu<-cluster[u]

    for(vv in g$edges[[uu]])
    {
      alt = dists[u] + d[uu,vv]
      v<-which(cluster==vv)
      if(alt < dists[v])
      {
        dists[v]<-alt
        previous[v]<-u
      }
      if(vv == j) left<-left=="B"
    }
    left[u]<-FALSE
    if(dbg)cat(paste("left: ",sum(left),"\r"))
  }
  if(dbg)cat("\n")

  path<-NULL
  jj<-which(j==cluster)
  u<-jj
  while(!is.na(previous[u]))
  {
    path<-c(path,cluster[u])
    u = previous[u]

  }
  path<-(rev(c(path,i)))
  path_length<-function(path,d){S<-0; for(i in 2:length(path)) S<- S + d[path[i-1],path[i]];return(S)}



  return(list(d=path_length(path,d), path=path))
}

#EOF


