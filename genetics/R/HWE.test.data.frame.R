
HWE.test.data.frame <- function(x, ..., do.Allele.Freq=TRUE,
                                do.HWE.test=TRUE)
{

  data <- makeGenotypes(x)
  names <- names(data)[sapply(data, is.genotype)]
  
  for(i in names)
  {
    gene   <- getlocus(i)
    genedata <- data[[i]]
    
    cat("\n")
    cat("+-------------------------------------\n");
    if(!is.null(gene))
      {
        cat("|\tMarker:\t ")
        print(gene)
      }
    else
      cat("|\tMarker: ", i, "\n")
    cat("+-------------------------------------\n");
  
    if(do.Allele.Freq)
      {
         # compute and print the allele and genotype frequencies
        sum  <-  summary(genedata)
        print(sum)
      }
  
    if(do.HWE.test)
      {
        if(length(allele.names(genedata))<2)
          {
            cat( '*** No variant alleles observed, unable to perform\n',
                 '*** test for Hardy-Wienburg Equilibrium. \n', sep='')
          }
        else
          {
            # now do and print the HWE test
            hwe  <- HWE.test(genedata, ...)
            print(hwe)
          }
  
      }
  }
  
}
