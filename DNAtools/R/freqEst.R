freqEst <- function(x){
  idcol <- seq(from=2,by=2,len=(ncol(x)-1)/2)
  x1 <- do.call("cbind",lapply(x[,idcol],paste))
  x2 <- do.call("cbind",lapply(x[,idcol+1],paste))
  names(x2) <- names(x1)
  xx <- as.data.frame(rbind(x1,x2))
  lapply(xx,function(allele) table(allele)/length(allele))
}
