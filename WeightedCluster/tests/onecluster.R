library(WeightedCluster)
data(mvad)
library(RUnit)

## Aggregating state sequence
aggMvad <- wcAggregateCases(mvad[1:10, 17:86], weights=mvad$weight[1:10])

## Creating state sequence object
mvad.seq <- seqdef(mvad[aggMvad$aggIndex, 17:86], weights=aggMvad$aggWeights)
## Computing Hamming distance between sequence
diss <- seqdist(mvad.seq, method="HAM")
clust4 <- wcKMedoids(diss, k=4, weights=aggMvad$aggWeights, cluster.only=TRUE)
dissd <- seqdist(mvad.seq, method="HAM", full.matrix=FALSE)
clust4d <- wcKMedoids(dissd, k=4, weights=aggMvad$aggWeights, cluster.only=TRUE)
#clust4d <- wcKMedoids(diss, k=4, weights=aggMvad$aggWeights,debuglevel=11)
checkEquals(clust4, clust4d)




## Compute the silhouette of each observation
clust1 <- rep(1, nrow(mvad.seq))
checkException(qual <- wcClusterQuality(diss, clust1, weights=aggMvad$aggWeights))

## Compute the silhouette of each observation
checkException(sil <- wcSilhouetteObs(diss, clust1, weights=aggMvad$aggWeights))

qual1 <- wcClusterQuality(diss, clust4, weights=aggMvad$aggWeights)
qual2 <- wcClusterQuality(dissd, clust4, weights=aggMvad$aggWeights)
checkEquals(qual1, qual2)

sil1 <- wcSilhouetteObs(diss, clust4, weights=aggMvad$aggWeights)
sil2 <- wcSilhouetteObs(dissd, clust4, weights=aggMvad$aggWeights)
checkEquals(sil1, sil2)
