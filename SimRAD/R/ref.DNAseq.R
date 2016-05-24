ref.DNAseq <-
function(FASTA.file, subselect.contigs = TRUE, prop.contigs= 0.10){
  ref <- readFasta(FASTA.file)
  ref <- sread(ref)
  N <- length(ref)
  if(N == 1){
    ref <- as.character(ref)
    return(ref)
  }
  if(N > 1 & subselect.contigs == FALSE){
    ref <- paste(ref, collapse="")
    return(ref)
  }
  if(N > 1 & subselect.contigs == TRUE){
    ref <- paste(c(ref[sample(1:length(ref), size=round(N*prop.contigs), replace=FALSE)]), collapse="")  
    return(ref)
  }
}
