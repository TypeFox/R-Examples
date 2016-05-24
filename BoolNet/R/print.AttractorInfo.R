# Custom print function for class AttractorInfo
print.AttractorInfo <- function(x, activeOnly = FALSE, ...)
{
  numGenes <- length(x$stateInfo$genes)
  attractors <- x$attractors
  
  for (i in seq_along(attractors))
  {
    if (is.null(attractors[[i]]$initialStates))
    # simple attractor
    {
      printSynchronousAttractor(getAttractorSequence(x, i), i, 
                                attractors[[i]]$basinSize, activeOnly=activeOnly)
    }
    else
    {
         # print general information on the attractor
      cat("Attractor ",i," is a complex/loose attractor consisting of ",ncol(attractors[[i]]$involvedStates),
        " state(s) and ",ncol(attractors[[i]]$initialStates), " transition(s)",sep="")
      
      if (activeOnly)
      {
        cat(".\nActive genes in the state transitions: \n")
        initialStates <- t(apply(attractors[[i]]$initialStates,2,function(state)
                               dec2bin(state,numGenes)))
        nextStates <- t(apply(attractors[[i]]$nextStates,2,function(state)
                               dec2bin(state,numGenes)))
                               
        binMatrix <- data.frame(initialStates,nextStates)                                 

        apply(binMatrix,1,function(row)
        {
          state1 <- paste(x$stateInfo$genes[which(row[seq_len(numGenes)] == 1)],collapse=", ")
          
          if (state1 == "")
            state1 <- "--"
          
          state2 <- paste(x$stateInfo$genes[which(row[seq_len(numGenes) + numGenes] == 1)],collapse=", ")
          
          if (state2 == "")
            state2 <- "--"
          
          cat(state1," => ",state2,"\n",sep="")
        })

      }
      else
      { 
        cat(":\n\n")
        initialStates <- apply(attractors[[i]]$initialStates,2,function(state)
                               paste(dec2bin(state,numGenes),collapse=""))
        nextStates <- apply(attractors[[i]]$nextStates,2,function(state)
                               paste(dec2bin(state,numGenes),collapse=""))

        binMatrix <- data.frame(initialStates,nextStates)

        apply(binMatrix,1,function(row)
        {         
            cat(row[1]," => ",row[2],"\n",sep="")
        })
        cat("\nGenes are encoded in the following order: ",paste(x$stateInfo$genes,collapse=" "),"\n\n",sep="")
      }      
    }
  }  
  return(invisible(x))
}
