ClusterPurity <-
function(ComputedClusters, TrueClasses) {
  sum(apply(table(ComputedClusters,TrueClasses), 1, max)) / length(ComputedClusters)
}
