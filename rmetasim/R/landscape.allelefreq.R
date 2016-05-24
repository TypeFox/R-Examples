#
# much faster allele freq calculator AES 3/6/09
#

landscape.allelefreq <- function (Rland, tbl.out = FALSE) 
{
    rv <- NULL
    for (i in 1:length(Rland$loci))
      {
        alleles <- landscape.locus(Rland,i)[, c(-1:-(landscape.democol()))]
#        print(i)
        if (landscape.ploidy(Rland)[i]==2)
          {
            adf <- as.data.frame(table(rep(landscape.populations(Rland),2),
                                       c(alleles[,1],alleles[,2])))
          } else {
            adf <- as.data.frame(table(landscape.populations(Rland),alleles))
          }
#        print(i)
        names(adf)[1:2] <- c("pop","alleles")
        adf$loc <- i
        popcnt <- aggregate(adf$Freq,by=list(pop=adf$pop),sum)
        adf <- merge(adf,popcnt)
        adf$Freq <- adf$Freq/adf$x #convert counts to proportions
        adf <- adf[,-5]
        adf <- adf[,c(2,3,1,4)] #just reorder the columns
        adf <- adf[adf$Freq>0,] #only include alleles present
        adf[order(adf$pop,adf$alleles),]
        rv <- rbind(rv,adf)
    }
    rownames(rv) <- 1:dim(rv)[1]
    if (tbl.out == TRUE) {
      xtabs(Freq ~ pop + alleles + loc, rv)
    }
    else {
      rv
    }
}
