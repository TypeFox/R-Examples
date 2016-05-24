# clustering has number of elements euqal the sum of numbers of rows in X and knowns
clusterAssigment <- function(X, knowns, clustering, method="euclidean") {
  # compute distances between X and knowns
  m = nrow(knowns)
  n = nrow(X)
  ppdist = as.matrix(dist(rbind(X, knowns)))[(n+1):(n+m),1:n]
  
  k = max(clustering)
  kkdist = matrix(NA, k, k, dimnames=list(1:k,1:k))
  
  clusteringX = clustering[1:n]
  clusteringk = clustering[(n+1):(n+m)]
  
  for (i in unique(clusteringX)) {
    for (j in unique(clusteringk)) {
      kkdist[j,i] = mean(ppdist[clusteringk==j, clusteringX==i, drop=F])
    } 
  } 

  # for every set of examples find closest cluster
  li = matrix(0,k,2)
  nk = length(unique(clusteringk))
  for (i in 1:nk) {
    px = which.min( apply(kkdist, 1, min) )
    py = which.min( kkdist[px, ] )
    li[i,] = as.numeric(c(names(px), names(py)))
    kkdist = kkdist[-px, -py, drop=F]
  }
  if (nk<k) {
    li[(nk+1):k,1] = setdiff(1:k, li[1:nk,1])
    li[(nk+1):k,2] = setdiff(1:k, li[1:nk,2])
  }
  #set new classes
  for (i in 1:k) {
    clustering[which(clusteringX==li[i,2])] = li[i,1]
  }  
  list(clustering = clustering, li = li)
}


# clustering has number of elements euqal the sum of numbers of rows in X and knowns
clusterAssigmentMeans <- function(knowns, class, means, method="euclidean") {
  # compute distances between X and knowns
  m = nrow(knowns)
  n = nrow(means)
  ppdist = as.matrix(dist(rbind(means, knowns)))[(n+1):(n+m),1:n]
  
  k = length(unique(class))
  kkdist = matrix(NA, k, nrow(means), dimnames=list(1:k,1:nrow(means)))
  for (i in 1:k) 
    kkdist[i, ] = colMeans(ppdist[class==i,, drop=F])

  # for every set of examples find closest cluster
  li = matrix(0,nrow(means),2)
  for (i in 1:k) {
    px = which.min( apply(kkdist, 1, min) )
    py = which.min( kkdist[px, ,drop=F] )
    li[i,] = as.numeric(c(colnames(kkdist)[py], rownames(kkdist)[px]))
    kkdist = kkdist[-px, -py, drop=F]
  }
  if (k < nrow(means)) {
     re = (k+1):nrow(means)
     li[re,2] = re
     li[re,1] = as.numeric(colnames(kkdist))
  }
  #set new classes
  list(li = li)
}


