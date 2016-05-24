#These are functions that manipulate the landscape to produce summary statistics
#implemented in the 'ape' or 'pegas' packages on CRAN.
#pegas and ape must be installed and loaded for these to work.
#interface to theta.h
landscape.theta.h <- function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            alleledist <- as.factor(landscape.locus(rland.tmp,lnum=j)[,c(-1:-(landscape.democol()))])
            if (length(unique(alleledist))>1)
              retval[i,j] <- pegas::theta.h(alleledist)
            else
              retval[i,j] <- NA
          }
      }
    retval
  }

#interface to theta.k
landscape.theta.k <- function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            alleledist <- as.factor(landscape.locus(rland.tmp,lnum=j)[,c(-1:-(landscape.democol()))])
            if (length(unique(alleledist))>1)
              retval[i,j] <- pegas::theta.k(alleledist)
            else
              retval[i,j] <- NA
          }
      }
    retval
  }

#waterson's segregating sites
landscape.theta.s<- function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            retval[i,j] <- NA
            if ((rland$loci[[j]]$type==253))
              {
                statevec <- landscape.locus.states(rland.tmp,lnum=j)$state
                seqlen <- nchar(statevec[1])
                                        #            print(j)
                print(paste("len statevec",length(statevec)))
                                        #            print(seqlen)
                segsites <- sum(apply(matrix(unlist(strsplit(statevec,split='')),ncol=seqlen,byrow=TRUE),2,function (x){length(unique(x))})>1)
            print(segsites)
                if (segsites>0)
                  retval[i,j] <- pegas::theta.s(segsites,seqlen)[1]
              }
          }
      }
    retval
  }

landscape.tajima.d <- function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            retval[i,j] <- NA
            if ((rland$loci[[j]]$type==253))
              {
                alleledist <- as.factor(landscape.locus(rland.tmp,lnum=j)[,c(-1:-(landscape.democol()))])
                if (length(unique(alleledist))>1)
                  theta.ewens <- c(theta.k(alleledist),0)
                else
                  theta.ewens <- c(NA,NA)
                statevec <- landscape.locus.states(rland.tmp,lnum=j)$state
                seqlen <- nchar(statevec[1])
                                        #            print(j)
                                        #                print(paste("len statevec",length(statevec)))
                                        #            print(seqlen)
                segsites <- sum(apply(matrix(unlist(strsplit(statevec,split='')),ncol=seqlen,byrow=TRUE),2,function (x){length(unique(x))})>1)
                if (segsites>0)
                  retval[i,j] <- (theta.ewens[1]-theta.s(segsites,seqlen)[1])/sqrt(theta.ewens[2]+theta.s(segsites,seqlen,variance=TRUE)[2])
              }
          }
      }
    retval
  }
