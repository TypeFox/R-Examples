stepArchetypesRawData <- function(data,numArch,numRep=3,verbose=TRUE){
  
  mycall <- match.call()
  as <- list()
  for (i in 1:length(numArch)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypes"
    for (j in seq_len(numRep)) {
      if (verbose) 
       cat("\n*** numArch=", numArch[i], ", rep=", j, ":\n", sep = "")
       as[[i]][[j]] <- archetypes(data, k = numArch[i], 
                                  family = archetypesFamily("original",scalefn = no.scalefn,
                                                                            rescalefn = no.rescalefn))
    }
  }
  return(structure(as, class='stepArchetypes',call=mycall))
}
