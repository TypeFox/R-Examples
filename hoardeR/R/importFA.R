# This function reads in a fasta file and prepares the vector from it
  importFA <- function(file){
    res <- readLines(file)
  # Check if the Fasta file is alternating, one line label, the next line sequence
    sumAlternating <- sum(grepl(">",res)==c(TRUE,FALSE))
    if(sumAlternating!=length(res)) stop("Your Fasta file is malformed. Please ensure that name rows start with > and that names and sequence
                                         rows are alternating.")
    seq <- res[seq(2,length(res),2)]
    names(seq) <- res[seq(1,length(res)-1,2)]
    class(seq) <- "fa"
    seq
} 
