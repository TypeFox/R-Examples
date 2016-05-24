rearrange <- function(genoGroups,rows.gex,genoSamples){
  newOrder <- c()
  for(i in 1:length(rows.gex))
  {
    ## WARNING!!! IN CASE OF MULTIPLE MEASURES IN THE GENOTYPE, WE PICK ALWAYS THE FIRST OCCURENCE!!!!
    newOrder[i] <- which((rows.gex[i]==genoSamples)==TRUE)[1]
  }

  ifelse(is.matrix(genoGroups) , output <- genoGroups[newOrder,], output <- genoGroups[newOrder])
  ifelse(is.matrix(genoGroups) , rownames(output) <- genoSamples[newOrder], names(output) <- genoSamples[newOrder])
  
  output
}