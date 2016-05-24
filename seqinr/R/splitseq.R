splitseq <- function(seq, frame = 0, word = 3){
#
# Compute all start positions of words to be returned:
#
  starts <- seq(from = frame + 1, to = length(seq), by = word)
#
# Extract them all:
#
  res <- sapply(starts, function(x) c2s(seq[x:(x + word - 1)]))
#
# remove last one if uncorrect length:
#
  if(nchar(res[length(res)]) != word) res <- res[-length(res)]
  return(res)
}
