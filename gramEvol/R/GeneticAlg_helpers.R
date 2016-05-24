ga.mutation <- function(genome, mutationChance, genomeLen = length(genome), 
                        genomeMin, genomeMax, allowrepeat,
                        dempeningFactor = 1) {
  
  mut_genomeLoc = runif(genomeLen) < mutationChance # do mutation for some of variables by chance
  num_muts = sum(mut_genomeLoc)
  
  # OPTION 1
  # mutate to something random
  #mutation = genomeMin[mut_genomeLoc] +
  #    runif(num_muts)*(genomeMax[mut_genomeLoc]-genomeMin[mut_genomeLoc]);
  
  # OPTION 2
  # mutate around solution
  direction       = runif(num_muts) - 0.5 # [-0.5 -> 0.5]
  mutationRange   = genomeMax[mut_genomeLoc]-genomeMin[mut_genomeLoc]
  mutation = round(genome[mut_genomeLoc] +  direction*mutationRange*dempeningFactor)
  # check if it is in domain. if not, then take random
  bad_mutations = which( (mutation < genomeMin[mut_genomeLoc]) | (mutation > genomeMax[mut_genomeLoc]) )
  for (b in bad_mutations) {
    mutation[bad_mutations] = ga.rand.int(n=1, 
                                          genomeMin[mut_genomeLoc][b],
                                          genomeMax[mut_genomeLoc][b])
  }
  
  # apply mutation
  genome[mut_genomeLoc] = mutation;
  if (!allowrepeat) {
    genome = ga.unique.maker(genome, genomeMin, genomeMax)
  }
  
  return (list(newGenome = genome, numMutations = num_muts))
}

ga.new.chromosome <- function(genomeLen, genomeMin, genomeMax, allowrepeat) {
  chromosome = round(runif(genomeLen) * (genomeMax - genomeMin) + genomeMin)
  
  if (!allowrepeat) {
    chromosome = ga.unique.maker(chromosome, genomeMin, genomeMax)
  }
  
  return (chromosome)
}

ga.rand.int <- function(n, mins, maxs) {
  (mins - 1) + sample.int(maxs - mins + 1, n, replace=TRUE)
}

ga.unique.maker <-
  function(x, genomeMin, genomeMax) {
    while (TRUE) {
      dup = duplicated(x)
      if (!any(dup))
        break
      for (i in which(dup)) {
        x[i] = x[i] + 1
        if (x[i] > genomeMax[i]) {
          x[i] = genomeMin[i]
        }
      }
    }
    
    return (x)
  }
