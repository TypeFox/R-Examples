h.func <-
function(geno, locustype, r, s, alleles){
  ## this is the top-level function for calculate h for a single individual
  ## geno: multilocus genotype
  ## r: parental allele frequencies of species with h=1
  ## s: parental allele frequencies of species with h=0
  maxH<-optimize(f=like.h, interval=0:1, maximum=TRUE,
                 tol=.Machine$double.eps^.5, geno=geno, locustype=locustype,
                 r=r, s=s, alleles=alleles)
  a<-0.0000001
  b<-maxH$maximum
  if(s.wrapper(a, b, maxH$objective, geno, locustype, r, s, alleles)){
    lowerBound<-uniroot(f=support.limit,
                        geno=geno, locustype=locustype, r=r, s=s, alleles=alleles,
                        max.obj=maxH$objective,
                        interval=c(a,b))$root
  } else{
    lowerBound<-0 #a
  }

  a<-maxH$maximum
  b<-1-0.0000001
  if(s.wrapper(a, b, maxH$objective, geno, locustype, r, s, alleles)){
    upperBound<-uniroot(f=support.limit,
                        geno=geno, locustype=locustype, r=r, s=s, alleles=alleles,
                        max.obj=maxH$objective,
                        interval=c(a,b))$root
  } else{
    upperBound<-1 #b
  }
  c(lowerBound, maxH$maximum, upperBound)
}

