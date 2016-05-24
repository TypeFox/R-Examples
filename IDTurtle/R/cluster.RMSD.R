cluster.RMSD <-
function(RMSD,method="ward.D"){
  RMSD<-as.matrix(RMSD)
  RMSD<-as.dist(RMSD)
  RMSDcluster<-hclust(RMSD,method)
  plot(RMSDcluster)
  return(RMSDcluster)
}
