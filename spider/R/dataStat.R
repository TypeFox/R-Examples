dataStat <- function(sppVector, genVector, thresh = 5){
unSpp <- unique(sppVector)
unGen <- unique(genVector)
sppnum <- sapply(unSpp, function(x) length(which(sppVector %in% x)))
tab <- table(NULL)
tab[1:7] <- c(length(unGen), length(unSpp), min(sppnum), max(sppnum), median(sppnum), mean(sppnum), length(which(sppnum < thresh)))
names(tab) <- c("Genera", "Species", "Min", "Max", "Median", "Mean", "Thresh" )
round(tab, digits=0)
}
