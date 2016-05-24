## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("kohonen")) {
  kohonen.present <- FALSE
  cat("Package kohonen not available - some code may not run.\nInstall it by typing 'install.packages(\"kohonen\")'")
} else {
  kohonen.present <- TRUE
}
if (!require("cluster")) {
  cluster.present <- FALSE
  cat("Package cluster not available - some code may not run.\nInstall it by typing 'install.packages(\"cluster\")'")
} else {
  cluster.present <- TRUE
}
if (!require("mclust")) {
  mclust.present <- FALSE
  cat("Package mclust not available - some code may not run.\nInstall it by typing 'install.packages(\"mclust\")'")
} else {
  mclust.present <- TRUE
}

data(wines, package = "ChemometricsWithRData")
wines.sc <- scale(wines)

## Hierarchical clustering
set.seed(7)
subset <- sample(nrow(wines), 20)
wines.dist <- dist(wines.sc[subset,])
wines.hcsingle <-  hclust(wines.dist, method = "single")

plot(wines.hcsingle, labels = vintages[subset])
wines.hccomplete <- hclust(wines.dist, method = "complete")
plot(wines.hccomplete, labels = vintages[subset])

wines.cl.single <- cutree(wines.hcsingle, h = 3.3)
table(wines.cl.single, vintages[subset])

wines.dist <- dist(wines.sc)
wines.hcsingle <- hclust(wines.dist, method = "single")
table(vintages, cutree(wines.hcsingle, k = 3))

wines.hccomplete <- hclust(wines.dist, method = "complete")
table(vintages, cutree(wines.hccomplete, k = 3))

if (cluster.present) {
  wines.agness <- agnes(wines.dist, method = "single")
  wines.agnesa <- agnes(wines.dist, method = "average")
  wines.agnesc <- agnes(wines.dist, method = "complete")
  cbind(wines.agness$ac, wines.agnesa$ac, wines.agnesc$ac)
}

## Partitional clustering
set.seed(21)
wines.km <- kmeans(wines.sc, centers = 3)
wines.km

table(vintages, wines.km$cluster)

## set.seed(3)
## wines.km <- kmeans(wines.sc, centers = 3)
## best <- wines.km
## for (i in 1:100) {
##   tmp <- kmeans(wines.sc, centers = 3)
##   if (sum(tmp$withinss) < sum(best$withinss))
##     best <- tmp
## }
## best

wines.km <- kmeans(wines.sc, centers = 3, nstart = 100)
wines.km

wines.pam <- pam(wines.dist, k = 3)
wines.pam

plot(wines.pam, main = "Silhouette plot")

best.pam <- pam(wines.dist, k = 2)
for (i in 3:10) {
  tmp.pam <- pam(wines.dist, k = i)
  if (tmp.pam$silinfo$avg.width < best.pam$silinfo$avg.width)
    best.pam <- tmp.pam
}
best.pam$medoids

table(vintages, best.pam$clustering)

if (mclust.present) {
  ## Probabilistic clustering
  wines.BIC <- mclustBIC(wines.sc, modelNames = "VVV")
  plot(wines.BIC)
  
  ## Error: map is the wrong one... Pay attention to the order the
  ## packages are loaded, since map is masked by something else
  wines.mclust2 <- mclustModel(wines.sc, wines.BIC)
  wines.mclust3 <- mclustModel(wines.sc, wines.BIC, G = 3)
  
  coordProj(wines.sc, dimens = c(7, 13),
            parameters = wines.mclust2$parameters,
            z = wines.mclust2$z, what = "classification")
  title("2 clusters: classification")
  coordProj(wines.sc, dimens = c(7, 13),
            parameters = wines.mclust3$parameters,
            z = wines.mclust3$z, what = "classification")
  title("3 clusters: classification")
  coordProj(wines.sc, dimens = c(7, 13),
            parameters = wines.mclust2$parameters,
            z = wines.mclust2$z, what = "uncertainty")
  title("2 clusters: uncertainty")
  coordProj(wines.sc, dimens = c(7, 13),
            parameters = wines.mclust3$parameters,
            z = wines.mclust3$z, what = "uncertainty")
  title("3 clusters: uncertainty")
  
  wines.BIC <- mclustBIC(wines.sc)
  plot(wines.BIC, legendArgs = list(x = "bottom", ncol = 2))
}

if (kohonen.present) {
  ## Comparing clusterings
  X <- scale(wines)
  set.seed(7)
  som.wines <- som(X, grid = somgrid(6, 4, "hexagonal"))
  set.seed(17)
  som.wines2 <- som(X, grid = somgrid(6, 4, "hexagonal"))
  som.hc <- cutree(hclust(dist(som.wines$codes)), k = 3)
  som.hc2 <- cutree(hclust(dist(som.wines2$codes)), k = 3)
  plot(som.wines, "mapping", bgcol = terrain.colors(3)[som.hc],
       pch = as.integer(vintages), main = "Seed 7")
  plot(som.wines2, "mapping", bgcol = terrain.colors(3)[som.hc2],
       pch = as.integer(vintages), main = "Seed 17")
  
  
  som.clust <- som.hc[som.wines$unit.classif]
  som.clust2 <- som.hc2[som.wines2$unit.classif]
  AdjRkl(som.clust, som.clust2)
  
  AdjRkl(som.clust, som.clust)
  
  AdjRkl(som.clust, som.clust2)
}
