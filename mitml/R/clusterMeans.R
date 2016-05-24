clusterMeans <- function(x, cluster, adj=FALSE, group=NULL){
# calculate cluster means

  # get objects if names are given
  isname <- c(length(x)==1, length(cluster)==1, length(group)==1) &
    c(is.character(x), is.character(cluster), is.character(group))
  if(any(isname)){
    parent <- parent.frame()
    if(isname[1]) x <- eval(parse(text=x),parent)
    if(isname[2]) cluster <- eval(parse(text=cluster),parent)
    if(isname[3]) group <- eval(parse(text=group),parent)
  }

  # prepare group
  if(!is.null(group)) {
    if(is.character(group)) group <- as.factor(group)
    if(is.factor(group)) group <- as.integer(group)
    ngr <- length(unique(group))
  }

  # format cluster (and groups)
  if(!is.numeric(cluster)) cluster <- as.integer(cluster)
  if(!is.null(group)) cluster <- cluster + group/(ngr+1)
  cluster <- match(cluster, unique(cluster)) 


  n.obs <- rowsum(as.integer(!is.na(x)), cluster)
  gm <- rowsum(x, cluster, na.rm = T)/n.obs
  gm[is.nan(gm)] <- NA
  gm <- gm[cluster]
  if(adj){
    n.obs <- n.obs[cluster]
    ((n.obs * gm) - x)/(n.obs - 1)
  }else{
    gm
  }
}
