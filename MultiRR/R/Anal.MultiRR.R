Anal.MultiRR <-
function(x){
All <- list()
All$SeriesVCV <- list()
All$IDVCV <- list()
All$pwI <- list()
All$pwS <- list()
All$n.obs <- list()
All$Individuals <- rep(NA, ncol(x$DataFrames))
All$Series <- rep(NA, ncol(x$DataFrames))
for(i in 1:ncol(x$DataFrames)){
w <- x$DataFrames[,i]

a <- lapply(w, lmerAll)

All$SeriesVCV[[i]] <- t(as.matrix(sapply(a, "[[", 1)))
All$IDVCV[[i]] <- t(as.matrix(sapply(a, "[[", 2)))
All$pwI[[i]] <- as.matrix(sapply(a, "[[", 3))
All$pwS[[i]] <- as.matrix(sapply(a, "[[", 4))
All$RInt[[i]] <- as.matrix(sapply(a, "[[", 5))
All$RSlope[[i]] <- as.matrix(sapply(a, "[[", 6))
All$Rrn[[i]] <-as.matrix(sapply(a, "[[", 7))
All$Residuals[[i]] <-as.matrix(sapply(a, "[[", 8))
All$Intercept[[i]] <- as.matrix(sapply(a, "[[", 9))
All$Slope[[i]] <-as.matrix(sapply(a, "[[", 10))
All$SimVCVInd[[i]] <- x$VCVInd
All$SimVCVSeries[[i]] <- x$VCVSeries
All$SimResiduals[[i]] <- x$Residuals
All$SimInt[[i]] <- x$SimInt
All$SimSlope[[i]] <- x$SimSlope
All$n.obs[[i]] <- nrow(as.data.frame(w[[1]]))
All$Individuals[i] <- mean(x$Individuals[,i])
All$Series[i] <- mean(x$Series[,i])
}
All$n.sim <- x$n.sim
All$n.ss <- ncol(x$DataFrames)
All
}
