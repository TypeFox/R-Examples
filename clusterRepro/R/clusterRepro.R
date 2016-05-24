permuteCol <- function(x) {
  dd <- dim(x)
  n <- dd[2]
  p <- dd[1]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(x[order(mm)], p, n, byrow = FALSE)
}

IGP.clusterRepro <- function(Data, Centroids) {
  Result <- c()
  Class <- rep(NA, ncol(Data))
  names(Class) <- colnames(Data)
  Features <- rownames(Data)
  Features <- Features[which(Features %in% rownames(Centroids))]
  Data <- Data[Features,]
  Centroids <- Centroids[Features,]
  Correlation <- cor(Data, Centroids, use = "pairwise.complete.obs", method = "pearson")
  for (i in 1:nrow(Correlation)) {
    Class[i] <- which.max(Correlation[i,])
  }
  Correlation <- cor(Data, method = "pearson", use = "pairwise.complete.obs")
  diag(Correlation) <- -Inf
  Classification <- matrix(0, ncol(Centroids), ncol(Centroids))
  for (i in 1:nrow(Correlation)) {
    Classification[Class[i], Class[which.max(Correlation[i,])]] <- Classification[Class[i], Class[which.max(Correlation[i,])]] + 1
  }
  Result$Class <- Class
  Result$IGP <- diag(Classification)/rowSums(Classification)
  Result$Size <- rowSums(Classification)
  return(Result)
}

clusterRepro <- function(Centroids, New.data, Number.of.permutations) {
  Result <- c()
  Features <- rownames(New.data)
  Features <- Features[which(Features %in% rownames(Centroids))]
  New.data <- New.data[Features,]
  Centroids <- Centroids[Features,]
  Calculation <- IGP.clusterRepro(New.data, Centroids)
  Actual.IGP <- Calculation$IGP
  Actual.Size <- Calculation$Size
  IGP <- c()
  Size <- c()

  SVD <- svd(Centroids)
  Centroids.prime <- Centroids%*%SVD$v
  for (i in 1:Number.of.permutations) {
    New.centroids.prime <- permuteCol(Centroids.prime)
    New.centroids <- New.centroids.prime%*%t(SVD$v)
    rownames(New.centroids) <- rownames(Centroids)
    Calculation <- IGP.clusterRepro(New.data, New.centroids)
    IGP <- c(IGP, Calculation$IGP)
    Size <- c(Size, Calculation$Size)
  }
  Result$p.value <- c(NA, length(Actual.Size))
  Result$Number <- c(NA, length(Actual.Size))
  for (i in 1:length(Actual.Size)) {
    if (length(which(Size == Actual.Size[i])) == 0 || Actual.Size[i] == 0) {
      Result$p.value[i] <- NA;
      Result$Number[i] <- 0;
    } else {
      Same.size <- IGP[which(Size == Actual.Size[i])];
      Result$p.value[i] <- length(which(Same.size >= Actual.IGP[i]))/length(Same.size);
      Result$Number[i] <- length(Same.size);
    }
  }
  Result$Actual.IGP <- Actual.IGP
  Result$Actual.Size <- Actual.Size
  return(Result)
}
