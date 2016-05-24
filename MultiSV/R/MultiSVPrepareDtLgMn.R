##' @export
`PrepMultiDtLgMn.default` <- function (MultiData){
  SampleSize = (dim(MultiData$MultiSvRw)[2]-3)
  MultiDataScl <- subset(MultiData$MultiSvRw,
  				(rowSums(MultiData$MultiSvRw[,4:dim(MultiData$MultiSvRw)[2]]))>=SampleSize & 
  				((rowSums(MultiData$MultiSvRw[,4:dim(MultiData$MultiSvRw)[2]])) <= SampleSize*1000) )
  colnames(MultiDataScl)[1] <- "chr";colnames(MultiDataScl)[2] <- "start";colnames(MultiDataScl)[3] <- "end"
  LocalScaleFactor <- as.character(unique(MultiDataScl[,1]))
  CrFc <- as.data.frame(matrix(nrow = length(LocalScaleFactor), ncol = dim(MultiDataScl)[2] -3) )
  row.names(CrFc) <- LocalScaleFactor
  colnames(CrFc) <- colnames(MultiDataScl[,4:dim(MultiDataScl)[2]])
  for(local in LocalScaleFactor){
    SubScale <- subset(MultiDataScl, MultiDataScl$chr==local)
    CrFc[local,] <- log(colMeans(SubScale[,4:(dim(MultiDataScl)[2])]))

  }  
  MultiData$MultiSvLgMn <- CrFc;rm(CrFc)
  MultiData
}




