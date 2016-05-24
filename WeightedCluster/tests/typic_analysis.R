library(WeightedCluster)
data(mvad)


## Aggregating state sequence
aggMvad <- wcAggregateCases(mvad[, 17:86], weights=mvad$weight)

## Creating state sequence object
mvad.seq <- seqdef(mvad[aggMvad$aggIndex, 17:86], weights=aggMvad$aggWeights)
## Computing Hamming distance between sequence
diss <- seqdist(mvad.seq, method="HAM")


averageClust <- hclust(as.dist(diss), method="average", members=aggMvad$aggWeights)

avgClustQual <- as.clustrange(averageClust, diss, weights=aggMvad$aggWeights, ncluster=20)

wardClust <- hclust(as.dist(diss), method="ward", members=aggMvad$aggWeights)

wardClustQual <- as.clustrange(wardClust, diss, weights=aggMvad$aggWeights, ncluster=20)

kmedClustQual <- wcKMedRange(diss, kvals=2:20, weights=aggMvad$aggWeights)

kmedwardClustQual <- wcKMedRange(diss, kvals=2:20, weights=aggMvad$aggWeights, initialclust=wardClust)


print(summary(avgClustQual, max.rank=3), digits=2)
print(summary(wardClustQual, max.rank=3), digits=2)
print(summary(kmedClustQual, max.rank=3), digits=2)
print(summary(kmedwardClustQual, max.rank=3), digits=2)

plot(wardClustQual, stat = c("ASWw", "HG", "PBC", "HC"))
avg4Medoids <- disscenter(diss, group=avgClustQual$clustering$cluster4, weights=aggMvad$aggWeights, medoids.index="first")

final <- wcKMedoids(diss, k=4, weights=aggMvad$aggWeights, initialclust=avgClustQual$clustering$cluster4)
final2 <- wcKMedoids(diss, k=4, weights=aggMvad$aggWeights, initialclust=avg4Medoids)
all.equal(final, final2)
final2


## Compute the silhouette of each observation
qual <- wcClusterQuality(diss, avgClustQual$clustering$cluster4, weights=aggMvad$aggWeights)

## Compute the silhouette of each observation
sil <- wcSilhouetteObs(diss, avgClustQual$clustering$cluster4, weights=aggMvad$aggWeights)


