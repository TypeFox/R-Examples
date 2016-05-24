 sil.score <- function (mat, nb.clus = c(2:13), nb.run = 100, iter.max = 1000, method = "euclidean") {

  sil.mean <- numeric()

  for (i in nb.clus) {

    avg.width <- numeric()

    j <- 1

    while (j <= nb.run) {

      cluster <- Kmeans(mat, i, iter.max = iter.max, method = method)$cluster

      len <- length(unique(cluster))

      #Kmeans may not return k specified clusters
      if (i == len) {

	#get silhouette for the same method
	silhouette <- silhouette(cluster, dist(mat, method))

	#get total mean of silhouette widths
	avg.width[j] <- as.numeric(summary(silhouette)[4])

	j <- j + 1
      }
    }
    #compute mean for each nb.clus
    sil.mean[i] <- mean(avg.width)
  }

return(sil.mean)
}

