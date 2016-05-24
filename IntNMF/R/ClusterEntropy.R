ClusterEntropy <-
function(ComputedClusters, TrueClasses) {
cross.tab <- table(ComputedClusters,TrueClasses)
  # Compute the inner sum for each clusters
inner.sum <- apply(cross.tab,1,function(x){
c.size <- sum(x)
sum(x * ifelse(x!=0,log2(x/c.size),0))})
# Calculate entropy
- sum(inner.sum)/(sum(cross.tab) * log2(ncol(cross.tab)))
}
