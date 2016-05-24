### R code from vignette source 'COMMUNAL.Rnw'

###################################################
### code chunk number 1: COMMUNAL.Rnw:41-52
###################################################
library(COMMUNAL)
data(BRCA.100)
varRange <- seq(20,100,20)
ks <- 2:8
measures <- c("average.between", "dunn", "widestgap", "dunn2", 
                      "pearsongamma", "g3", "max.diameter", "avg.silwidth")

BRCA.results <- clusterRange(dataMtx=BRCA.100, ks = ks,
                              varRange=varRange, 
                              validation=measures,
                              verbose = F)


###################################################
### code chunk number 2: COMMUNAL.Rnw:58-61
###################################################
algs <- getGoodAlgs(BRCA.results, algs="all")

algs


###################################################
### code chunk number 3: COMMUNAL.Rnw:66-73
###################################################
monotoneClusterRange(BRCA.results)

measuresCorr(BRCA.results)

measures <- getNonCorrNonMonoMeasures(BRCA.results, goodAlgs=algs, numMeasures = 4)

measures


###################################################
### code chunk number 4: COMMUNAL.Rnw:79-83
###################################################
plot.data <- plotRange3D(BRCA.results, ks, algs, measures, plot3D=F)

## print the values from the 3D plot
plot.data


###################################################
### code chunk number 5: COMMUNAL.Rnw:96-99
###################################################
result <- BRCA.results$all.results$vars_20
clusters <- result$getClustering(k=3)
apply(clusters, 2, table)


###################################################
### code chunk number 6: COMMUNAL.Rnw:103-111
###################################################
# re-key cluster labels to most frequent assignments
mat.key <- clusterKeys(clusters)
examineCounts(mat.key)

# find 'core' clusters
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters
head(core) # the cluster assignments


###################################################
### code chunk number 7: COMMUNAL.Rnw:115-120
###################################################
clusters.example <- data.frame(
  alg1=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,1)),
  alg2=as.integer(c(1,1,1,1,1,3,3,3,3,3,2,2,2,2,1)),
  alg3=as.integer(c(3,3,3,3,3,1,1,1,1,1,2,2,2,2,2))
)


###################################################
### code chunk number 8: COMMUNAL.Rnw:124-126
###################################################
mat.key <- clusterKeys(clusters.example)
mat.key # cluster indices are relabeled


###################################################
### code chunk number 9: COMMUNAL.Rnw:130-131
###################################################
examineCounts(mat.key)


###################################################
### code chunk number 10: COMMUNAL.Rnw:135-137
###################################################
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters


###################################################
### code chunk number 11: COMMUNAL.Rnw:141-143
###################################################
core <- returnCore(mat.key, agreement.thresh=99)
table(core) # 0 is undetermined


