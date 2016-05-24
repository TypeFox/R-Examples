insilico.digest <- function (DNAseq, cut_site_5prime1, cut_site_3prime1, cut_site_5prime2="NULL",  cut_site_3prime2="NULL",  
                               cut_site_5prime3="NULL", cut_site_3prime3="NULL", cut_site_5prime4="NULL", cut_site_3prime4="NULL", verbose = TRUE) 
{
    recognition_code1 <- paste(cut_site_5prime1, cut_site_3prime1, sep = "")
    digest1 <- insilico.digest.internal(DNAseq, recognition_code1, cut_site_5prime1, cut_site_3prime1)
    if(cut_site_5prime2 == "NULL"){RESULT <- digest1}
    if(cut_site_5prime2 != "NULL"){
       recognition_code2 <- paste(cut_site_5prime2, cut_site_3prime2, sep = "")
       digest2 <- insilico.digest.internal(digest1, recognition_code2, cut_site_5prime2, cut_site_3prime2)
       RESULT <- digest2
    }
    if(cut_site_5prime3 != "NULL"){
       recognition_code3 <- paste(cut_site_5prime3, cut_site_3prime3, sep = "")
       digest3 <- insilico.digest.internal(digest2, recognition_code3, cut_site_5prime3, cut_site_3prime3)
       RESULT <- digest3
    }
    if(cut_site_5prime4 != "NULL"){
       recognition_code4 <- paste(cut_site_5prime4, cut_site_3prime4, sep = "")
       digest4 <- insilico.digest.internal(digest3, recognition_code4, cut_site_5prime4, cut_site_3prime4)
       RESULT <- digest4
    }
    RESULT <- unlist(RESULT)
    if(verbose == TRUE){
      if(cut_site_5prime2 == "NULL"){
        cat("Number of restriction sites: ", length(insilico.digest.internal(DNAseq, recognition_code1, cut_site_5prime1, cut_site_3prime1))-1, "\n", sep="")
      } 
      if(cut_site_5prime2 != "NULL" && cut_site_5prime3 == "NULL"){
      cat("Number of restriction sites for the first enzyme: ", length(insilico.digest.internal(DNAseq, recognition_code1, cut_site_5prime1, cut_site_3prime1))-1, "\n", sep="")
      cat("Number of restriction sites for the second enzyme: ", length(insilico.digest.internal(DNAseq, recognition_code2, cut_site_5prime2, cut_site_3prime2))-1, "\n", sep="")
      dig1 <- RESULT[isMatchingStartingAt(cut_site_3prime1, RESULT)]
      dg1 <- reverseComplement(DNAStringSet(dig1))
      re2match <- reverseComplement(DNAStringSet(cut_site_5prime2))
      dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
      dig2 <- reverseComplement(dg2)
      dig1bis <- RESULT[isMatchingStartingAt(cut_site_3prime2, RESULT)]
      dg1bis <- reverseComplement(DNAStringSet(dig1bis))
      re2matchbis <- reverseComplement(DNAStringSet(cut_site_5prime1))
      dg2bis <- dg1bis[isMatchingStartingAt(re2matchbis[[1]], dg1bis)]
      dig2bis <- reverseComplement(dg2bis)
      cat("Number of type AB and BA fragments:", length(dig2)+length(dig2bis), "\n", sep="")
      RE1RE1.re2match <- reverseComplement(DNAStringSet(cut_site_5prime1))
      RE1RE1.dg2 <- dg1[isMatchingStartingAt(RE1RE1.re2match[[1]], dg1)]
      RE1RE1.dig2 <- reverseComplement(RE1RE1.dg2)
      cat("Number of type AA fragments:", length(RE1RE1.dig2), "\n", sep="")
      dig3 <- RESULT[isMatchingStartingAt(cut_site_3prime2, RESULT)]
      dg3 <- reverseComplement(DNAStringSet(dig3))
      RE2RE2.re2match <- reverseComplement(DNAStringSet(cut_site_5prime2))
      RE2RE2.dg3 <- dg3[isMatchingStartingAt(RE2RE2.re2match[[1]], dg3)]
      RE2RE2.dig3 <- reverseComplement(RE2RE2.dg3)
      cat("Number of type BB fragments:", length(RE2RE2.dig3), "\n", sep="")
      }
      if(cut_site_5prime3 != "NULL"){
        cat("Number of restriction sites 1: ", length(insilico.digest.internal(DNAseq, recognition_code1, cut_site_5prime1, cut_site_3prime1))-1, "\n", sep="")
        cat("Number of restriction sites 2: ", length(insilico.digest.internal(DNAseq, recognition_code2, cut_site_5prime2, cut_site_3prime2))-1, "\n", sep="")
        cat("Number of restriction sites 3: ", length(insilico.digest.internal(DNAseq, recognition_code3, cut_site_5prime3, cut_site_3prime3))-1, "\n", sep="")        
      }
      if(cut_site_5prime4 != "NULL"){
        cat("Number of restriction sites 1: ", length(insilico.digest.internal(DNAseq, recognition_code1, cut_site_5prime1, cut_site_3prime1))-1, "\n", sep="")
        cat("Number of restriction sites 2: ", length(insilico.digest.internal(DNAseq, recognition_code2, cut_site_5prime2, cut_site_3prime2))-1, "\n", sep="")
        cat("Number of restriction sites 3: ", length(insilico.digest.internal(DNAseq, recognition_code3, cut_site_5prime3, cut_site_3prime3))-1, "\n", sep="") 
        cat("Number of restriction sites 4: ", length(insilico.digest.internal(DNAseq, recognition_code4, cut_site_5prime4, cut_site_3prime4))-1, "\n", sep="")               
      }            
    } 
    return(RESULT)
}