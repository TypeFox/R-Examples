translate <- function(seq, frame = 0, sens = "F", numcode = 1, NAstring = "X",
ambiguous = FALSE)
{

  if(any(seq%in%LETTERS)){
    seq <- tolower(seq)
  }
  #
  # Take the reverse complementary strand when required:
  #

  if(sens == "R") seq <- comp(rev(seq), ambiguous = ambiguous)

  #
  # Transform the sequence in its numerical encoding equivalent
  # with textbook order, that is t = 0, c = 1, a = 2, g = 3
  #

  seqn <- s2n(seq, levels = s2c("tcag"))

  #
  # Compute the length of the sequence when its length in codons
  # is an integer:
  #
  
  l <- 3*((length(seq) - frame) %/% 3)

  #
  # Compute the indices for the first codon positions:
  #
  
  c1 <- seq(from = frame + 1, to = frame + l, by = 3)
  
  #
  # Compute the indices of codons in the translation table:
  #

  tra <-  16*seqn[c1] + 4*seqn[c1 + 1] + seqn[c1 + 2] + 1

  #
  # Get the translation table:
  #

  code <- s2c(SEQINR.UTIL$CODES.NCBI$CODES[numcode])

  #
  # Translate the sequence:
  #

  result <- code[tra]

  #
  # Replace missing values by the string for missing amino-acids:
  #

  result[is.na(result)] <- NAstring
  
  #
  # More work is required if ambiguous bases are handled:
  #
  
  if(ambiguous){
    toCheck <- which(result == NAstring)
    for( i in toCheck ){
      codon <- seq[c1[i]:(c1[i]+2)]
      allcodons <- as.vector(outer(as.vector(outer(amb(codon[1]), amb(codon[2]), paste, sep = "")), amb(codon[3]), paste, sep = ""))
      allaminoacids <- sapply(allcodons, function(x) translate(s2c(x), numcode = numcode, ambiguous = FALSE))
      if( all(allaminoacids == allaminoacids[1])) result[i] <- allaminoacids[1]
    }
  }

  return( result )
}

