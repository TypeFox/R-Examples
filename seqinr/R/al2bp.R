al2bp <- function(allele.name, repeat.bp = 4, offLadderChars = "><", split = "\\."){
  allele.name <- as.character(allele.name)
  if(any(s2c(allele.name) %in% s2c(offLadderChars))) return(NA)	
  dec <- unlist(strsplit(allele.name, split = split))
  res <- repeat.bp*suppressWarnings(as.numeric(dec[1]))
         # NA are returned for X and Y at Amelogenin locus for instance
  if(length(dec) > 1) res <- res + as.numeric(dec[2])
  return(res)
}
