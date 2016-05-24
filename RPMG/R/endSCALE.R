endSCALE<-function(arange, digits = 3)
  {
   ## print(arange)

    arange = range(arange)

    ### adjust the digits if there are not enough to handle the difference
    dd =diff(arange)
    fdd = format(dd, scientific = TRUE)
    fss = strsplit(fdd, "e") 
     Dexp = abs(  as.numeric( sapply(fss, "[[", 2) ) )

    if( digits<Dexp)
      {
        digits=Dexp+1

      }
    #######  
    char = format(arange,digits = digits ,  scientific = TRUE)
    ss = strsplit(char, "e") 
    Aexp = as.numeric( sapply(ss, "[[", 2) )
    Mant = as.numeric( sapply(ss, "[[", 1) )
    MAexp = max(Aexp)
    NewMantdigits = Mant*10^(Aexp-MAexp)

    ###   here need to check for very small values
    ##  newchar = paste(NewMant, "*x*10^{", MAexp , "}", sep = "" )
    
    newexp = paste("x*10^{", MAexp , "}", sep = "" )
    return(c( as.character(NewMantdigits ) , newexp  ) )
  }
