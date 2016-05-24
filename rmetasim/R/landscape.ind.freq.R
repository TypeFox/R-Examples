
landscape.ind.freq <- function(Rland,include.states=TRUE)
  {
      l <- Rland
      aml <- vector("list",length(landscape.ploidy(l)))
      for (loc in 1:length(landscape.ploidy(l)))
      {
          genos <- landscape.locus(l,loc)[,-1:-landscape.democol()]
          ploidy <- landscape.ploidy(l)[loc]
          if (l$loci[[loc]]$type!=253)
          {
              lst <- landscape.locus.states(l,loc)
              names(lst$state) <- lst$aindex
              if (ploidy==2)
              {
                  genos[,1] <- unname(lst$state[as.character(genos[,1])])
                  genos[,2] <- unname(lst$state[as.character(genos[,2])])
              } else {
                  genos <- unname(lst$state[as.character(genos)])
              }
          }
          amat <- sapply(as.numeric(names(table(genos))),function(x,genos,pl)
          {
              if (pl==2)
              {
                  (genos[,1]==x)+(genos[,2]==x)
              } else
              {
                  genos==x
              }
          },genos=genos,pl=ploidy)
          aml[[loc]] <- apply(amat,2,function(x,pl){x/pl},pl=ploidy) #allele freqs per ind
      }
      do.call(cbind,aml)
  }
