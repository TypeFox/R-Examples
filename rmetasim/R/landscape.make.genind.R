#Create genind object from landscape
#ignores the haploid loci
landscape.make.genind <- function(Rland)
  {
      tab <- landscape.ind.freq(Rland)*2
      dimnames(tab) <- list(rownames=1:dim(tab)[1],colnames=landscape.freq.locnames(Rland))
      pl <- landscape.ploidy(Rland)
      populations <- landscape.populations(Rland)
      gi=genind(tab,pop=as.factor(populations),ploidy=2)
      gi[,loc=which(pl>1)]
  }
