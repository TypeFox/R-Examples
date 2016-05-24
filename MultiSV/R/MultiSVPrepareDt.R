##' @export
`PrepMultiDt.default` <- function (MultiData){
  SampleSize = (dim(MultiData$MultiSvRw)[2]-3)
  MultiDataScl <- subset(MultiData$MultiSvRw,
  (rowSums(MultiData$MultiSvRw[,4:dim(MultiData$MultiSvRw)[2]]))>=SampleSize & ((rowSums(MultiData$MultiSvRw[,4:dim(MultiData$MultiSvRw)[2]])) <= SampleSize*1000) )
  colnames(MultiDataScl)[1] <- "chr";colnames(MultiDataScl)[2] <- "start";colnames(MultiDataScl)[3] <- "end"
  LocalScaleFactor <- as.character(unique(MultiDataScl[,1]))
  CrFc <- as.data.frame(matrix(nrow = 1, ncol = dim(MultiDataScl)[2]))
  MltScl <- MultiDataScl
  for(local in LocalScaleFactor){
    SubScale <- subset(MultiDataScl, MultiDataScl$chr==local)
      for (i in 4:(dim(MultiDataScl)[2])){
    		CrFc[1,i] <- (mean(colSums(SubScale[,4:(dim(MultiDataScl)[2])])))/(sum(SubScale[,i]))
    		}
    	for (i in 4:(dim(MultiDataScl)[2])){
    		MltScl[rownames(SubScale),i] <- SubScale[,i] * CrFc[1,i]
    	}
  }
  MultiData$MultiSvScl <- MltScl;rm(MltScl);rm(MultiDataScl)
  MultiData
}

