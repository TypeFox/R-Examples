######Transform grid into adj matrix
Grid.Adjmatrix.Transfer=function(grid, euclidean.dist=1)
{
  net=dist(grid,diag=TRUE, upper=TRUE)
  net=as.matrix(net)
  net=sapply(1:nrow(net), function(i,j) return(net[i,j]<=euclidean.dist))
  net=net*1
  diag(net)<-0
  return(net)
}
