Fis <- function(tab.pop, allele.column)
  {
    
    # Calculates allele or genotype frequencies from tab.pop NA 0 etc are not omitted
    # Inputformat
    # Individual   population allele.column allele.column+1 ....
    #   x              EGB         1             2
    #   y              EGB         2             4
    #   z              EGB         3             5
    #   .               .          .             .
  
    # NA's are excluded for each locus
    tab.pop <- data.frame(tab.pop[,1],tab.pop[,2],tab.pop[,allele.column],tab.pop[,allele.column+1])
    tab.pop <- tab.pop[complete.cases(tab.pop),]	
    popsize <- length(tab.pop[,1])

    # Transform from alleles to genotypes
    
    tab.freq.gen <- data.frame(tab.pop[,1],tab.pop[,2],paste(tab.pop[,3],tab.pop[,4],sep="-"))    
    tab.freq.gen <- table(tab.freq.gen[,3])/sum(table(tab.freq.gen[,3]))
    
    sum.homo <- sum(tab.freq.gen[which(lapply(lapply(strsplit(as.character(names(tab.freq.gen)), "-", fixed = TRUE),as.numeric),diff)==0)])
    
    # Building and aggregating tables
       
    tab.list <- list(table(tab.pop[,3]), table(tab.pop[,4]))
    
    tab.freq <- do.call(rbind, lapply(lapply(tab.list, unlist), "[", unique(unlist(c(sapply(tab.list,names))))))
    colnames(tab.freq) <- unique(unlist(c(sapply(tab.list,names))))


    # Calculating frequencies
        
    tab.freq <- colSums(tab.freq, na.rm=TRUE)/sum(tab.freq, na.rm=TRUE)
    tab.freq <- tab.freq[which(names(tab.freq)!=0)]
    
    if (sum.homo==1){fis <- 1}
    
    if (sum.homo!=1)
    {
    Hsb <- popsize/(popsize-1)*(1-sum(tab.freq*tab.freq)-((1-sum.homo)/2/popsize))
    fis <- 1-((1-sum.homo)/(Hsb))
    }
    
    tab.total.weir <- weir(tab.pop,tab.freq,popsize)
    
    # Arithmetic mean of alleles for Weir and Cockerham Fis
    # Empirical weighted Weir Fis for each locus

    b.weir <- sum(tab.total.weir[2,])
    c.weir <- sum(tab.total.weir[3,])

    weir.locus <- 1-(sum(c.weir)/(sum(c.weir+b.weir)))
    
    if (is.nan(weir.locus)==TRUE){weir.locus <- 1}

    fis.weir <- list(tab.total.weir,weir.locus)
    names(fis.weir) <- c("Allele Information","Weighted mean of allele Fis")
    
    tab.freq[order(names(tab.freq))]
    tab.freq.gen[order(names(tab.freq.gen))]
    
    fis.return <- list(tab.freq, tab.freq.gen, sum.homo, popsize, fis, fis.weir)
    names(fis.return) <- c("Frequency_of_Alleles","Frequency_of_Genotypes","Homozygosity", "Population_size", "Fis_Nei", "Fis_Weir")
    return(fis.return)
    
  } 

