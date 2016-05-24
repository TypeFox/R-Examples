write.fasta <- function(sequences, names, file.out, open = "w", nbchar = 60, 
    as.string = FALSE){
  #
  # Open output file:
  #
  outfile <- file(description = file.out, open = open)
  
  #
  # Function to write one sequence in output file:
  #
  write.oneseq<-function(sequence, name, nbchar, as.string){
    writeLines(paste(">", name, sep = ""), outfile)
    if(as.string) sequence <- s2c(sequence)
    l <- length(sequence)
    q <- floor(l/nbchar)
    r <- l - nbchar*q
    if(q > 0){
      sapply(seq_len(q), function(x) writeLines(c2s(sequence[(nbchar*(x - 1) + 1):(nbchar*x)]), outfile))
    }
    if(r > 0){
      writeLines(c2s(sequence[(nbchar*q + 1):l]), outfile)
    }
  }
  
  #
  # Write all sequences in output file:
  #
  if(!is.list(sequences)){
    write.oneseq(sequence = sequences, name = names, nbchar = nbchar, as.string = as.string)
  } else {
    n.seq <- length(sequences)
    sapply(seq_len(n.seq), function(x) write.oneseq(sequence = as.character(sequences[[x]]), 
      name = names[x], nbchar = nbchar, as.string = as.string))
  }
  #
  # Close output file:
  #
  close(outfile)
}
