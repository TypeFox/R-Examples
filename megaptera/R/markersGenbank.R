markersGenbank <- function(taxa, classify = FALSE){
  
  id <- esearch(taxa, db = "Nucleotide")
  x <- lapply(id, EFetchLocus)
  if ( classify ){
    classifyLoci <- function(x){
      x <- gsub("large subunit ribosomal RNA|ribosomal RNA|rRNA", "rDNA", x)
      x <- gsub("translation elongation factor 1 alpha", "tef1", x)
      if ( length(grep("contains", x)) > 0 ){
        x <- unlist(strsplit(x, ", | and "))
        x <- gsub("^contains |^and ", "", x)
      }
      x <- gsub("rDNA gene", "rDNA", x)
      x <- gsub("internal transcribed spacer ", "ITS", x)
      x
    }
    x <- sapply(x, classifyLoci)
  }
  names(x) <- id
  x
}