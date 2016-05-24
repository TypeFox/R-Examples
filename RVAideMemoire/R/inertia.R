inertia <- function(tab,w=rep(1,nrow(tab)),center=TRUE,adj=FALSE) {
  n <- nrow(tab)
  m <- if (center) {
    apply(tab,2,function(x) {wmean(x,w=w)})
  } else {
    rep(0,ncol(tab))
  }
  ef <- ifelse(adj,n-1,n)
  res <- sum(apply(tab,1,function(x) {sum((x-m)^2)}))/ef
  return(res)
}
