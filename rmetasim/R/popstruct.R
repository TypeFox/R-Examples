#
#Population structure stats
#

#samples the entire landscape, does not correct for sample bias
#on the assumption that every individual is sampled

landscape.Fst <- function(rland,verb=FALSE)
  {
    aft <- landscape.allelefreq(rland,TRUE)
    Fst <- matrix(0,length(landscape.ploidy(rland)),length(aft[1,,1]))
    for (i in 1:length(Fst[,1]))
      {
        Fst[i,] <- ((colSums(sweep(aft[,,i],2,colMeans(aft[,,i]))^2)/(rland$intparam$habitats))/(colMeans(aft[,,i])*(1-colMeans(aft[,,i]))))
      }
    if (verb)
      {
        print(paste("Populations:",rland$intparam$habitats,"Loci:",rland$intparam$locusnum))
        print(paste("Mean per locus:"))
        print(rowMeans(Fst,na.rm=TRUE))
        print(paste("Overall mean:",mean(rowMeans(Fst,na.rm=TRUE))))
        print("===")
      }
    Fst
  }
