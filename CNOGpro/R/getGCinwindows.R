getGCinwindows <-
function(sequence, windowlength){
  gcvector <- numeric()
  q <- 1
  condition <- T
  chromosome <- sequence[[1]]
  
  while(condition){
    gcvector <- c(gcvector, GC(chromosome[q:(q+(windowlength-1))]))
    q <- q + windowlength
    if((length(chromosome)-q) <= windowlength){
      condition <- F
    }
  }
    
  gcvector <- c(gcvector, GC(chromosome[q:length(chromosome)]))
  return(gcvector)
}
