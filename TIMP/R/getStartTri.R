"getStartTri" <- function(xx, nexp, ncomp){
  y <- as.matrix(rowSums(xx))
  dim(y) <- c(ncomp, nexp)
  ampList <- vector("list", length=nexp)
  for(i in 1:nexp) ampList[[i]] <- rep(0,ncomp)
  fixedamps <- vector("list", length=ncomp) 
  for(i in 1:ncomp) {
    ycol <- y[i,]
    mind <- which(ycol==max(ycol))
    ycol <- ycol /  max(ycol)
    fixedamps[[i]] <- c(mind,i)
    for(j in 1:nexp) 
      ampList[[j]][i] <- ycol[j]
  }
  list(fixedamps=fixedamps, ampList=ampList)
}
