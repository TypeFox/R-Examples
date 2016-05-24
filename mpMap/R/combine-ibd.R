combine_ibd <- function(obsgeno)
{
  geno <- list()
  nmrk <- ncol(obsgeno$founders)/2

  dom <- apply(obsgeno$founders[,1:nmrk], 2, 
	 function(x) return(sum(names(table(x))==c("0","1"))==2))
  snp <- apply(obsgeno$founders[,1:nmrk], 2, 
	 function(x) return(sum(names(table(x))==c("0","2"))==2))

  geno$founders <- obsgeno$founders[,1:nmrk]
#  geno$founders[,dom] <- geno$founders[,dom]+obsgeno$founders[,(nmrk+which(dom==TRUE))]
#  geno$founders[,dom][geno$founders[,dom]>0] <- 1
#  geno$founders[,snp] <- (geno$founders[,snp]+obsgeno$founders[,(nmrk+which(snp==TRUE))])/2

  geno$finals <- obsgeno$finals[,1:nmrk]
  geno$finals[obsgeno$finals[,1:nmrk]!=obsgeno$finals[,nmrk+1:nmrk]] <- NA
#  geno$finals[,dom] <- geno$finals[,dom]+obsgeno$finals[,(nmrk+which(dom==TRUE))]
#  geno$finals[,dom][geno$finals[,dom]>0] <- 1
#  geno$finals[,snp] <- (geno$finals[,snp]+obsgeno$finals[,(nmrk+which(snp==TRUE))])/2

  return(geno)
}
