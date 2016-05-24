# LAST UPDATE: 2013-11-21

revCompTest <- function(DNAbin, split, benchmark = TRUE,
                        mafft.exe){
  
  if ( benchmark ) if ( !split["benchmark"] ) split <- !split
  DNAbin <- as.list(DNAbin)
  rc <- function(s) rev(as.DNAbin(comp(s)))
  DNAbin[!split] <- lapply(DNAbin[!split], rc)
  mafft(DNAbin, path = mafft.exe)
}
