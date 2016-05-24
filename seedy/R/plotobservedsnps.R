plotobservedsnps <- function(data, timepoint=1, coverage=50, error=0.001, iterations=100, maxsnp=50, legend=TRUE, 
                    ylim=c(0,1.5*coverage), ...) {
  
  if (!"obs.freq"%in%names(data)) {
    stop("Data object must be simulated under 'full' sampling")
  }
  
  genos <- data$obs.strain[[timepoint]]
  genos.fr <- data$obs.freq[[timepoint]]
  str <- data$librstrains
  snplist <- NULL
  snpnum <- NULL
  for (i in 1:length(genos)) {
    geno_num <- which(str == genos[i]) #this is the key to extracting SNP info from the libr object
    muts <- data$libr[[geno_num]]
    if (sum(is.na(muts))==0) {
      #nucs <- libr$nucs[[geno_num]]
      for (j in 1:length(muts)) {
        #tag <- paste(nucs[1],muts[1],sep = "-") - removing b/c it is not important to account for multiple mutations at the same site
        tagpos <- snplist == muts[j]
        if (sum(tagpos) > 0) {      # there must be a nicer way than this!
          tagnum <- which(tagpos)
          snpnum[tagnum] = snpnum[tagnum]+genos.fr[i] ## IMPORTANT - need to check to see if this holds under all growth conditions
        }
        else {      
          snplist <- append(snplist,muts[j])
          snpnum <- append(snpnum,genos.fr[i])
        }
      }
    }
  }
  
  eq.size <- sum(genos.fr)
  snpnum[snpnum > eq.size] <- eq.size # this fixes the multiple mutation issue.  See notes at top.

  #######
  
  gcov <- rpois(eq.size,coverage)
  polys <- rbinom(eq.size, gcov, error)
  
  #######
  
  polys <- polys[polys > 0]
  polys <- polys[order(polys, decreasing = TRUE)]
  
  ######
  
  snp.sim <-data.frame(NULL)
  snp.freq <- snpnum/eq.size  
  for (i in 1:iterations){
    res.vec <-rep(0,length(snpnum))
    for (j in 1:length(snpnum)){      
      #number of true SNP covered
      t.snps <- rpois(1,(snp.freq[j]*coverage))
      # number of true non-SNPs
      t.nonsnps <- rpois(1,((1-snp.freq[j])*coverage))
      #false postivie
      false.pos <-rpois(1,t.nonsnps*error)
      #number of real SNPs that err back to the true base (ie er/3)
      false.neg <-rpois(1,t.snps*error/3)
      res.vec[j] <- t.snps+false.pos-false.neg  
    }
    snp.sim <-rbind(snp.sim,res.vec)
  }
  colnames(snp.sim) <- snplist

  #######
  
  snp.sim <- snp.sim[order(colMeans(snp.sim), decreasing = TRUE)]
  
  if (maxsnp > length(colMeans(snp.sim))) {
    maxsnp <- length(colMeans(snp.sim))
  }
  if (maxsnp > length(polys)) {
    maxsnp <- length(polys)
  }
  
  plot(NULL, xlim=c(1, maxsnp), ylim=ylim, ...)
  
  for (j in 1:maxsnp) {
    lines(c(j,j), c(quantile(snp.sim[,j], c(0.025,0.975))))
    points(j,mean(snp.sim[,j]), pch=16, col="blue")
  }
  lines(polys[1:maxsnp], col = "red")
  abline(h = coverage, lty = 3)
  if (legend) {
    legend("topright",c("Neutral SNPs","Sequence error","Av. coverage"), col = c("blue","red","black"),
           pch = c(16,NA,NA), lty = c(1,1,3))
  }
  return(invisible(snp.sim))
}