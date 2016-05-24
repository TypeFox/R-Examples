adapt.select <-
function(sequences, type="AB+BA", cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, cut_site_3prime2){
  if(type == "AA"){
    dig1 <- sequences[isMatchingStartingAt(cut_site_3prime1, sequences)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(cut_site_5prime1))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    return(dig2)
  }
  if(type == "BB"){
    dig1 <- sequences[isMatchingStartingAt(cut_site_3prime2, sequences)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(cut_site_5prime2))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    return(dig2)    
  }
  if(type == "AB"){
    dig1 <- sequences[isMatchingStartingAt(cut_site_3prime1, sequences)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(cut_site_5prime2))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    return(dig2)    
  }   
  if(type == "BA"){ 
    dig1 <- sequences[isMatchingStartingAt(cut_site_3prime2, sequences)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(cut_site_5prime1))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    return(dig2)      
  }
  if(type == "AB"){
    dig1 <- sequences[isMatchingStartingAt(cut_site_3prime1, sequences)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(cut_site_5prime2))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    return(dig2)    
  }
  if(type == "AB+BA"){
    dig1 <- sequences[isMatchingStartingAt(cut_site_3prime1, sequences)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(cut_site_5prime2))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    dig1bis <- sequences[isMatchingStartingAt(cut_site_3prime2, sequences)]
    dg1bis <- reverseComplement(DNAStringSet(dig1bis))
    re2matchbis <- reverseComplement(DNAStringSet(cut_site_5prime1))
    dg2bis <- dg1bis[isMatchingStartingAt(re2matchbis[[1]], dg1bis)]
    dig2bis <- reverseComplement(dg2bis)
    dig3 <- c(dig2, dig2bis)        
    return(dig3)      
  } 
}
