thin<-function(dat,win=100,shift=25,r2thresh=0.01) { 
  # write .ped and .map files, call plink to do thinning. plink returns
  # a list of rsnumbers of snps to keep. Subset full list of snps accordingly
  pedinf<-unrel.pre.pedinf(nrow(dat$snp.data))
  write.pedfile(pedinf,dat$snp.data,file="mydata.ped")
  snpIDs<-colnames(dat$snp.data)
  write.mapfile(dat$snp.support,snpIDs,file="mydata.map")
  plink.comm<-paste("plink --noweb --file mydata --indep-pairwise",
              win, shift, r2thresh,sep=" ")
  system(plink.comm)
  prune.in<-scan("plink.prune.in",what=character())
  unlink(c("mydata.ped","mydata.map","plink.log","plink.prune.in",
           "plink.prune.out"))
  keep<-(snpIDs %in% prune.in)
  dat$snp.data<-dat$snp.data[,keep]
  dat$snp.support<-dat$snp.support[keep,]
  return(dat)
}

write.mapfile<-function(snp.support,snpIDs,file="out.map") {
     # from plink documentation, need a file with columns
     # chromosome (1-22, X, Y or 0 if unplaced)
     # rs# or snp identifier
     # Genetic distance (morgans)
     # Base-pair position (bp units)
     out<-data.frame(snp.support$Chromosome,snpIDs,
                     snp.support$Gen_loc*100, #Gen_loc is in cM
                     snp.support$Position)
    write.table(out,file,quote=FALSE,row.names=FALSE,col.names=FALSE) 
}     

unrel.pre.pedinf<-function (n) 
{
    ones <- rep(1, n)
    zeros <- rep(0, n)
    ids <- (1:n)
    pedinf <- cbind(ids, ids, zeros, zeros, ones, ones)
    return(pedinf)
}

