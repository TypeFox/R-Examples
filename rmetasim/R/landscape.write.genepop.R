#
#
# A quick R function to take a landscape and output a genepop file
# AES 12/1/2015
#
# This function will just ignore haploid loci and only output diploids.  It also uses the three digit allele designations
# It takes all individuals and puts them into the output file
#
# IMPORTANT to avoid the missing data indicator in GENEPOP, 1 is added to all allele numbers
#
# l is a landscape and fn is the name of the output file
#
landscape.write.genepop <- function(Rland,fn="genepop.out",title="rmetasim landscape output")
  {
    if (is.landscape(Rland)) #input error check
      {
        cat(paste(date(),": ",title,sep=''),file=fn)
        cat("\n",file=fn,append=T)
        #first output the locus names
        for (i in which(landscape.ploidy(Rland)==2)) #only use the diploids
          {
            cat(paste("locus-",i," ",sep=''),file=fn,append=T)
            cat("\n",file=fn,append=T)
          }
        #now outputting the genotypic data is a bit more complex
        #first make a big matrix of the genotypes is GENEPOP format
        #the rows still correspond to the rows in l$individuals
        printmat <- matrix("",ncol=(1+length(which(landscape.ploidy(Rland)==2))+1),
                           nrow=length(landscape.populations(Rland)))
        printmat[,1] <- sprintf("Class-%d Gen-%d,",Rland$individuals[,1],Rland$individuals[,3])
        colcount <- 2
        for (i in which(landscape.ploidy(Rland)==2))
          {
            if (Rland$loci[[i]]$type!=253)
              {
                printmat[,colcount] <- sprintf("%03d%03d ",landscape.states(Rland,i)[,(landscape.democol()+1)]+1,
                                               landscape.states(Rland,i)[,(landscape.democol()+2)]+1)
              }
            else
              {
                printmat[,colcount] <- sprintf("%03d%03d ",landscape.locus(Rland,i)[,(landscape.democol()+1)]+1,
                                               landscape.locus(Rland,i)[,(landscape.democol()+2)]+1)
              }
            colcount <- colcount + 1
          }
        printmat[,colcount] <- rep("\n",length(landscape.populations(Rland)))
        #now output one population at a time
        for (i in unique(landscape.populations(Rland)))
          {
            cat("POP\n",file=fn,append=T)
            cat(t(printmat[landscape.populations(Rland)==i,]),file=fn,append=T)
          }
      } else { #see, it really wasn't a landscape
        print ("make sure that the landscape is in the right form")
      }
  }
