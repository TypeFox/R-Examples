# COMMUNAL package
# part 3: extract core clusters

#### Create the test matrix model ####
flipLevels <- function(x, levels) {
  lev1<-levels(x)[levels[1]]
  lev2<-levels(x)[levels[2]]
  x[which(x==lev1)] <- NA
  x[which(x==lev2)] <- lev1
  x[which(is.na(x))] <- lev2
  x
}

induceError <- function(x, percent) {
  index <- sample(1:length(x), size=ceiling(length(x)*percent))
  x[index] <- sample(levels(x), size=length(index), replace=T)
  x
}

genClustMtx <- function(x, cols, errorPct=0.2) {
  mat <- x
  for(i in 1:(cols-1)){
    temp <- x
    for(i in 1:cols){ temp <- flipLevels(temp, as.numeric(sample(levels(temp), size=2))) }
    temp <- induceError(temp, errorPct)
    mat <- cbind(mat,temp)
  }
  mat
}

genTestMtx <- function(k, algs, rows, errorPerAlg=0.2) {
  x<-as.factor(sample(1:k, size=rows, replace=T))
  mat <- genClustMtx(x, cols=algs, errorPct=errorPerAlg)
  rownames(mat) <- seq(1:rows)
  colnames(mat) <- paste0("alg.", seq(1:algs))
  mat
}

#### Find duplicate rows, then find keys, then sort by keys ####
splitOutDuplicates <- function(mat, k) {
  clust <- list()
  i=1
  while (TRUE) {
    ident <- which(duplicated(mat))
    matches <- apply(mat, 1, function(x) identical(as.integer(x), as.integer(mat[ident[1],])))
    eval(parse(text=paste0("clust$cl",i," <- mat[matches, ,drop=FALSE]")))
    mat <- mat[!matches, ,drop=FALSE]
    if (length(setdiff(ident, which(matches))) == 0) break;
    i=i+1    
  }
  clust$remaining <- mat
  clust
}

getKeys <- function(uniqueClustList, k) {
  # Fails if less than k-1 unique groupings in the data
  lengths <- lapply(uniqueClustList, length)
  index <- order(unlist(lengths[-length(lengths)]), decreasing=TRUE)
  keys <- t(data.frame(uniqueClustList[[ index[1] ]][1,]))
  i=2
  while (TRUE) {
    clust.i.row <- uniqueClustList[[ index[i] ]][1,]
    if(!(clust.i.row[1] %in% keys[,1])) keys <- rbind(keys, clust.i.row);
    if(i==length(index) | dim(keys)[1]==k) break
    i<-i+1
  }
  if (dim(keys)[1]==(k-1)){ keys <- genRow(keys,k) }
  keys
}


swapKeys <- function(x, swapMtx) {
  x <- swapMtx[match(x, swapMtx[,2]), 1]
  x
}


genRow <- function(keys,k) {
  missing <- apply(keys, 2, function(x) setdiff(seq(1:k), x)[1])
  keys <- rbind(keys, unlist(missing))
  keys
}


reorderMat <- function(mat, index.alg) {
  if(!(index.alg %in% colnames(mat))) stop("index.alg not found in column names")
  x <- which(colnames(mat) == index.alg)
  reindex.mat <- cbind(mat[ , x, drop=F], mat[ ,-x, drop=F] )
  reindex.mat
}


clusterKeys.internal <- function(mat, k, index.alg=colnames(mat)[1]) {
  stopifnot(typeof(mat)=='integer')
  mat <- reorderMat(mat, index.alg)
  clust <- splitOutDuplicates(mat, k)
  if(length(clust)<(k-1)) stop("Not enough unique rows to robustly de-key matrix.\n Need more samples, fewer algorithms, or fewer clusters.\n")
  keys <- getKeys(clust,k)
  if(dim(keys)[1]<(k-1)){stop("The only unique keys available were:\n", apply(keys,1, function(x) c(x,"\n")))}
  # Next re-key according to unique rows
  if (dim(keys)[2] > 1) {
    for(j in 2:dim(keys)[2]){
      swap <- keys[ ,c(1,j)]  
      col <- lapply(mat[,j], function(x) swap[match(x, swap[,2]), 1])
      mat[,j] <- unlist(col)
    }  
  }
  mat[is.na(mat)] <-NA
  mat
}


optAlg <- function(clusters, k) {
  tmp <- vector()
  for(i in 1:length(colnames(clusters))){
    mat.key <- clusterKeys.internal(clusters, k, index.alg=colnames(clusters)[i])
    tmp <- c(tmp, examineCounts(mat.key)[1,2])
  }
  opt.Alg <- colnames(clusters)[match(min(tmp), tmp)]
  opt.Alg
}

cleanClusters <- function(clusters, k, min.size=3) {
  # remove poor clusterings (cluster too small or wrong number of clusters)
  clusters <- data.frame(clusters)
  goodnum <- lapply(lapply(clusters, table), length) == k
  if (!all(goodnum)) {
     cat(paste0(paste(names(which(!goodnum)), sep=", "), ": returned < ", k, " clusters - algorithm(s) dropped\n"))
  }
  goodsize <- lapply(lapply(clusters, table), min) >= min.size
  if (!all(goodsize)) {
    cat(paste0(paste(names(which(!goodsize)), sep=", "), ": returned cluster with size < ", min.size, " - algorithm(s) dropped\n"))
  }
  goodalgs <- goodsize & goodnum
  clean.clusters <- clusters[, goodalgs, drop=F]
  data.matrix(clean.clusters)
}

#### Public functions ####
clusterKeys <- function(clusters, min.size=3) {
  k <- max(clusters, na.rm=T)
  clusters <- cleanClusters(clusters, k, min.size)
  opt.Alg <- optAlg(clusters, k)
  mat.key <- clusterKeys.internal(clusters, k, opt.Alg)
  mat.key
}

examineCounts <- function(mat.key) {
  tmp <- table(apply(mat.key, 1, function(row) max(table(row)/length(row))))
  tmp <- cbind(percent.agreement=100*(as.numeric(names(tmp))), sample.counts=tmp)
  rownames(tmp) <- NULL
  percent.remaining <- vector()
  for(i in 1:dim(tmp)[1]) {
    percent.remaining <- c(percent.remaining, signif((sum(tmp[,2])-sum(tmp[1:i, 2]))/sum(tmp[,2]) ,2)) 
  }
  tmp <- cbind(tmp, percent.remaining.if.removed=percent.remaining)
  tmp
}

# v1.1 changed <= agreement.thresh to < agreement.thresh
returnCore <- function(mat.key, agreement.thresh=50) {
  if(agreement.thresh < 50) warning("agreement threshold <50; any ties resolved arbitrarily\n")
  # Here we examine whether the highest-voted cluster has less agreement than threshold
  # If less agreement, set to 0 
  index <- apply(mat.key, 1, function(row) max(table(row)/length(row)) <= agreement.thresh/100 )
  thresholded <- mat.key
  thresholded[index, ] <- 0
  core <- apply(thresholded, 1, function(row) {
    tab <- table(row)
    names(tab)[tab==max(tab)][1]
  })
  cat("A total of", sum(index), 'samples were rejected as not robustly clustered.\n')
  cat(signif(sum(!index)/length(index),3)*100, "% samples remain.\n")
  core
}
