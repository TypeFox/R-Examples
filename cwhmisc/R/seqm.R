seqm <- function(from, to, by=1) {
  if ((to-from)*by < 0) return(NULL)
  else return(seq(from, to, by))
}
