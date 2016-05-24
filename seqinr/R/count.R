count <- function(seq, wordsize, start = 0, by = 1, freq = FALSE, alphabet = s2c("acgt"), frame = start){
#
# For backward compatibility:
#
  if(!missing(frame)) start = frame
#
# istarts contains the first position of oligomers in the sequence (starting at 1)
#
  istarts <- seq(from = 1 + start, to = length(seq), by = by)
#
# oligos contains the first character of oligomers:
#
  oligos <- seq[istarts]
#
# oligos.levels contains all possible oligomers for a given alphabet:
#
  oligos.levels <- levels(as.factor(words(wordsize, alphabet = alphabet)))
#
# For n-mers with n >= 2 we paste the following characters in the
# sequence to build the observed set of all oligomers. Some NA are
# generated at the end of the sequence and discarded when counting
# them.
#
  if (wordsize >= 2){
    for(i in 2:wordsize){
      oligos <- paste(oligos, seq[istarts + i - 1], sep = "")
    }
  }
#
# We count all oligomers, even missing ones, and discard NA
#
  counts <- table(factor(oligos, levels = oligos.levels))
#
# Build result:
#
  if(freq == TRUE) counts <- counts/sum(counts)
  return(counts)
}
