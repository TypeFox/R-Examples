#Function: ghap.fst
#License: GPLv3 or later
#Modification date: 22 Jan 2016
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Calculate Fst each haploblock 

ghap.fst <- function(blockstats.pop1, blockstats.pop2, blockstats.tot){
  
  #Merge dataframes
  if(identical(blockstats.pop1$BLOCK,blockstats.pop2$BLOCK) & identical(blockstats.pop1$BLOCK,blockstats.tot$BLOCK)){
    freq <- NULL
    freq$BLOCK <- blockstats.tot$BLOCK
    freq$CHR <- blockstats.tot$CHR
    freq$BP1 <- blockstats.tot$BP1
    freq$BP2 <- blockstats.tot$BP2
    freq$EXP.H.pop1 <- blockstats.pop1$EXP.H
    freq$EXP.H.pop2 <- blockstats.pop2$EXP.H
    freq$EXP.H.tot <- blockstats.tot$EXP.H
    freq <- as.data.frame(freq)
  }else{
    stop("The three input dataframes should contain the same haplotype blocks\n")
  }
  
  #Compute FST
  Hs <- (freq$EXP.H.pop1+freq$EXP.H.pop2)/2
  freq$FST <- (freq$EXP.H.tot - Hs)/freq$EXP.H.tot
  freq$FST[which(freq$FST < 0)] <- 0
  freq$FST[which(freq$FST > 1)] <- 1
  
  #Return output
  return(freq)

}

