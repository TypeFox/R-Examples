### AES 1/22/06
### supposed to calculate pairwise relatedness among all individuals in a landscape
### reference: @ARTICLE{queller:1989a,
###  author = {Queller, D.C. and Goodnight, K.F.},
###  title = {Estimating relatedness using genetic markers},
###  journal = {Evolution},
###  year = {1989},
###  volume = {43},
###  pages = {258-275},
###}
###

landscape.relate <- function(rland)
  {
    ##calculate the allele frequecies for each locus ignoring the population designations
    ##haploid loci are not removed here, but will not be used in relatedness calcs
    acdf <- landscape.allelecount(rland)
    acnp <- aggregate(acdf$Freq,by=list(loc=acdf$loc,allele=acdf$allele),sum)
    acnp$tot <- rep(dim(rland$individuals)[1],dim(acnp)[1]) * landscape.ploidy(rland)[acnp$loc]
    acnp$Freq <- acnp$x/acnp$tot
    acnp$loc <- as.numeric(as.character(acnp$loc))
    acnp$allele <- as.numeric(as.character(acnp$allele))

    internalrelate(rland,acnp)
  }

internalrelate <- function(rland, acnp)
  {
    #indmat must be individuals locus data.  diploid only, with no classification data
    #acnp should have 4 columns: 1) locus 2) allele 3) count 4) total
    #ref individual is the 'x' ind in relate and part is the 'y' individual
    
    indmat <- as.matrix(rland$individuals[,-1:-(landscape.democol())])
    acnp <- as.matrix(acnp[,1:4])
    .Call("relateinternal",indmat,acnp,PACKAGE = "rmetasim")
  }

