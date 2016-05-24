like.h <-
function(geno, locustype, r, s, alleles, h){
  n.loci<-length(geno)
  result<-0
  ##  maxalleles<-dim(r)[2]
  for(i in 1:n.loci){
    ## two alleles at a locus, grab those, and parental allele
    ## frequencies and pass those along to per.locus.like
    locus.contents<-unlist(strsplit(as.character(geno[i]), "/"))
    allele1<-locus.contents[1]
    allele2<-locus.contents[2]
    if(locustype[i] == "C" || locustype[i] == "c"){
      a1.index<-which(alleles[i,] == allele1)
      a2.index<-which(alleles[i,] == allele2)
      ## make sure we have numeric indexes for both alleles, otherwise
      ## one of the alleles does not exist in the parents and does not
      ## contribute to the estimation of the hybrid index
      if(length(a1.index) == 0 || length(a2.index) == 0){
        next
      }
      else{
        result <- sum(result, per.locus.like(allele1, allele2, locustype[i],
                                          r[i, a1.index],
                                          r[i, a2.index],
                                          s[i, a1.index],
                                          s[i, a2.index], h),na.rm=TRUE)
      }
    }
    else if(locustype[i] == "D" || locustype[i] == "d"){
      a1.index<-which(alleles[i,] == allele1)
      ## need logic to get the index of the alternative allele, since
      ## there is just one in the genotype. a1.index should be 1 or 2
      ## at this point
      a2.index<-ifelse(a1.index==1, 2, 1)
      result <- sum(result, per.locus.like(allele1, allele2, locustype[i],
                                        r[i, a1.index],
                                        r[i, a2.index],
                                        s[i, a1.index],
                                        s[i, a2.index], h),na.rm=TRUE)
    }
    else { ## haploid locus
      a1.index<-which(alleles[i,] == allele1)
      result <- sum(result, per.locus.like(allele1, allele2, locustype[i],
                                        r[i, a1.index],
                                        r2=NULL,
                                        s[i, a1.index],
                                        s2=NULL, h),na.rm=TRUE)
    }
  }
  result
}

