sampleChromosome <-
function(chromosome){
  asample <- matrix(,nrow=20, ncol=20)
  for(i in 1:20){
    newsample <- sample(chromosome$ReadsprWindow, size=20, replace=T)
    asample[i,] <- newsample
  }
  chrmean <- mean(asample)
  # Must only sample in areas where CN=1. The following line of code assumes that samples with the higher variance are more
  # likely to include observations from duplications. Random samples from duplicated areas will greatly increase variance,
  # and we are trying to infer the variance in single-copy areas.
  chrvar <- median(apply(asample, MARGIN=1, var))
  
  return(list(chrmean, chrvar))
}
