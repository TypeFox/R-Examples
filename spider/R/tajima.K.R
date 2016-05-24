tajima.K <-
function(DNAbin, prop = TRUE){
   res <- mean(dist.dna(DNAbin, model="N"))
   if(prop) res <- res/dim(DNAbin)[2]
   res
}

