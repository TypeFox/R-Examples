siloF <-
function(oldSilo, Nbro, rank) {
  lengOld <- length(oldSilo)
  seqSilo <- seq(oldSilo[1], oldSilo[2], length = (Nbro + 1))
  out <- c(seqSilo[rank], seqSilo[rank + 1])
  return(out)
}
