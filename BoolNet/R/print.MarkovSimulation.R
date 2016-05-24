# Custom print function for class MarkovSimulation
print.MarkovSimulation <- function(x, activeOnly = FALSE, ...)
{
  genes <- x$genes
  numGenes <- length(genes)

  cat("States reached at the end of the simulation:\n")
  if (activeOnly)
  {
    reachedStates <- data.frame(apply(x$reachedStates[,seq_len(numGenes)],1,function(state)
                         {
                            r <- paste(genes[which(state == 1)],collapse=", ")
                            if (r == "")
                              r <- "--"
                            r
                         }),x$reachedStates$Probability)
    colnames(reachedStates) <- c("Active genes", "Probability")                        

    print(reachedStates)                           
  }
  else
  {
    print(x$reachedStates)
  }

  if (!is.null(x$table))
  {
                           
    cat("\nProbabilities of state transitions in the network:\n")
    transitionProbs <- getTransitionProbabilities(x)

    colIndices <- c(1,numGenes,numGenes + 1, 2*numGenes, 
                    2*numGenes + 1)
                      
    if(activeOnly)
    {
      inputStates <- apply(transitionProbs,1,function(row)
                          {
                            r <- paste(genes[which(row[colIndices[1]:colIndices[2]] == 1)],collapse=", ")
                            if (r == "")
                              r <- "--"
                            r
                          })
      outputStates <- apply(transitionProbs,1,function(row)
                            {
                              r <- paste(genes[which(row[colIndices[3]:colIndices[4]] == 1)],collapse=", ")
                              if (r == "")
                                r <- "--"
                              r
                            })
      colWidth <- max(c(sapply(inputStates,nchar),sapply(outputStates,nchar)))
      align <- "left"
    }
    else
    {
        inputStates <- apply(transitionProbs,1,function(row)
                          paste(row[colIndices[1]:colIndices[2]],collapse=""))
        outputStates <- apply(transitionProbs,1,function(row)
                          paste(row[colIndices[3]:colIndices[4]],collapse=""))  
        colWidth <- numGenes
        align <- "right"                  
    }

    binMatrix <- data.frame(inputStates,outputStates,
                            transitionProbs[,colIndices[5]])                    
                    
    cat(format("State",width=max(7,colWidth),justify=align),"    ",
        format("Next state",width=max(11,colWidth + 2),justify=align),
        format("Probability",width=13,justify="right"),"\n",sep="")                    
                    
    apply(binMatrix,1,function(row)
    {
      # paste all states of input and output into one string, and put out all columns of the table in a
      # formatted way
      cat(format(row[1],width=max(7,colWidth),justify=align),
          " => ",
          format(row[2],width=max(11,colWidth + 2),justify=align),
          format(row[3],width=13,justify="right"),"\n",sep="")
    })  
  }
  return(invisible(x))
}
