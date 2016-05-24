 kmeans.run <- function (mat, nb.clus = 2, nb.run = 1000, iter.max = 10000, method = "euclidean") {

  #result will be stored
  #in a list of two lists
  res <- list()

  res.kmeans <- list()

  #compute k-means according to nb.run
  j <- 1
  while (j <= nb.run) {
    cluster <- Kmeans(mat, nb.clus, iter.max = iter.max, method = method)
    len <- length(unique(cluster$cluster))
    #Kmeans function may not 
    #return k specified clusters
    if (nb.clus == len) {
      res.kmeans[[j]] <- cluster$cluster
      j <- j + 1
    }
  }
  
  #use ref to remove
  #redundant k-means results
  ref <- rep(1, nb.run)

  names(ref) <- c(1:nb.run)

  res.match <- list()

  for (i in 1:(nb.run - 1)) {
    #avoid useless calculations
    if (!is.na(match(i, names(ref)))) {

      res.match[[i]] <- list()

      for (j in (i + 1):nb.run) {
        #avoid useless calculations
        if (!is.na(match(j, names(ref)))) {

          tab <- table(res.kmeans[[i]], res.kmeans[[j]])
          #keep result of matchClasses for next calculations
          res.match[[i]][[j]] <- matchClasses(tab, verbose = FALSE, method = "exact")

          value <- sum(diag(tab[, res.match[[i]][[j]]])) / sum(tab)
          #if value is 1, it means that
          #the two k-means give the cluster attributions
          if (value == 1) {
            #count same results
            ref[match(i, names(ref))] <- ref[match(i, names(ref))] + 1
            #remove redundant results
            ref <- ref[-match(j, names(ref))]
          }
        }
      }
    }
  }

  #keep only maximum occurrence between k-means results
  index.max <- which.max(ref)
  value.max <- ref[index.max]
  ref1 <- ref[-index.max]

  k <- as.numeric(names(index.max))
  occ <- lapply(res.kmeans[[k]], function(i) {
    v <- rep(0, nb.clus); names(v) <- seq_len(nb.clus); v[i] <- value.max; return (v)})
  for (i in 1:length(ref1)) {
    l <- as.numeric(names(ref1[i]))
    if (k > l) {
      m <- res.match[[l]][[k]]
      mname<-as.vector(m)
      m <- as.numeric(names(m))
      names(m) <- mname
     } else {
      m <- res.match[[k]][[l]]
     }
    #match clus attribution between unique k-means results
    clus <- match(res.kmeans[[l]], m)
    for (j in 1:length(clus)) {
      #attribute occurrence to respective clus
      occ[[j]][clus[j]] <- occ[[j]][clus[j]] + ref1[i]
    }
  }


  occ.perc <- lapply(occ, function(i) {(i * 100) / nb.run})
  #keep only maximum occurrence between clus
  occ.perc.max <- lapply(occ.perc, function(i) {i[which.max(i)]})
  res$elements <- occ.perc

  #revert occ.perc.max
  n <- sapply(occ.perc.max, names)
  m <- sapply(occ.perc.max, c, use.names = FALSE)
  res$clusters <- tapply(m, n, c)
  class (res) <- c("kmean")
  return(res)
}

