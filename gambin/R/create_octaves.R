create_octaves <-
function(abundances, subsample = 0)  
{
	if(!subsample == 0) abundances <- .sample_abundances(abundances, subsample)
  stopifnot(is.numeric(abundances))
  abundances = abundances[abundances > 0]     # remove zeros
  octs = floor(sapply(abundances, log2))
  octs = factor(octs, levels = 0:max(octs))   # ensure that all octaves are tabled, even if no species fall in that octave
  ret <- data.frame(table(octs))
  names(ret) = c("octave","species")
  ret$octave = 0:(nrow(ret)-1)                # octaves are numbered from 0, 1, 2... etc.
  return(ret)
}
