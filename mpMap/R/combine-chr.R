combine_chr <- function(map, map.function)
{
  # combine chromosomes together into an object which has the probability 
  # of recombination between adjacent loci
  map2 <- unlist(map)
  
  # computes recombination fractions from map distances using cMfx
  genome <- vector()

  for (i in 1:length(map))
	genome <- c(genome, sapply(diff(map[[i]])/100, map.function), 0.5)

  return(genome[1:(length(genome)-1)])
}
