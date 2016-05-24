##' @export
`Bin2MultiSV.default` <- function(CfgFile){
  OutDt <- list()
  DataSt <- list()
  Conf <- read.table(CfgFile,header=F)
	colnames(Conf) <- c("Contrst","Id","Path")
	OutDt$Conf <- Conf
	for (i in 1:(dim(Conf)[1])){
		DataSt[[i]] <- read.table(sub("$", "", Conf[i,3]) ,header=T)
	}
	SortedTotal <- DataSt[[1]]
  for (i in 2:(length(DataSt))){
    SortedTotal <- merge(SortedTotal,(DataSt[[i]]),by=c("chr","start","end"),sort=T)
  }
	OutDt$MultiSvRw <- SortedTotal[order(SortedTotal$chr,SortedTotal$start),]
  OutDt
}