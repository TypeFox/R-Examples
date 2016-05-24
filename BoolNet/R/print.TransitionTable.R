# Custom print function for class TransitionTable
print.TransitionTable <- function(x, activeOnly=FALSE, ...)
{

  geneCols <- setdiff(colnames(x),c("attractorAssignment","transitionsToAttractor"))
  numGenes <- (length(geneCols)) / 2

  colIndices <- c(1,numGenes,numGenes + 1, 2*numGenes);
  
  if ("attractorAssignment" %in% colnames(x))
    colIndices <- c(colIndices, 2*numGenes + 1)
    
  if ("transitionsToAttractor" %in% colnames(x))
    colIndices <- c(colIndices, 2*numGenes + 2)    
      
  genes <- sapply(colnames(x)[seq_len(numGenes)],function(n)strsplit(n,".",fixed=TRUE)[[1]][2])
  
  if(activeOnly)
  {
    inputStates <- apply(x,1,function(row)
                        {
                          r <- paste(genes[which(row[colIndices[1]:colIndices[2]] == 1)],collapse=", ")
                          if (r == "")
                            r <- "--"
                          r
                        })
    outputStates <- apply(x,1,function(row)
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
      inputStates <- apply(x,1,function(row)
                        paste(row[colIndices[1]:colIndices[2]],collapse=""))
      outputStates <- apply(x,1,function(row)
                        paste(row[colIndices[3]:colIndices[4]],collapse=""))  
      colWidth <- numGenes
      align <- "right"                  
  }
  
  binMatrix <- cbind(inputStates,outputStates)
  
   if ("attractorAssignment" %in% colnames(x))
    binMatrix <- cbind(binMatrix, x[,colIndices[5]])
    
   if ("transitionsToAttractor" %in% colnames(x))
    binMatrix <- cbind(binMatrix, x[,colIndices[6]])

  binMatrix <- as.data.frame(binMatrix)

  cat(format("State",width=max(7,colWidth),justify=align),"    ",
      format("Next state",width=max(11,colWidth + 2),justify=align),
      if ("attractorAssignment" %in% colnames(x))
      {
        format("Attr. basin",width=13,justify="right")
      }
      else
        "",
      if ("transitionsToAttractor" %in% colnames(x))
      { 
        format("# trans. to attr.",width=19,justify="right")
      }
      else
        "",
      "\n",sep="")
      
  apply(binMatrix,1,function(row)
  {
    # paste all states of input and output into one string, and put out all columns of the table in a
    # formatted way
    cat(format(row[1],width=max(7,colWidth),justify=align),
        " => ",
        format(row[2],width=max(11,colWidth + 2),justify=align),
        if ("attractorAssignment" %in% colnames(x))
        { 
          format(row[3],width=13,justify="right")
        }
        else 
          "",
        if ("transitionsToAttractor" %in% colnames(x))
        { 
            format(row[4],width=19,justify="right")
        }
        else
        "",
        "\n",sep="")
  })  
  if (!activeOnly)
    cat("\nGenes are encoded in the following order: ",
        paste(genes,collapse=" "),"\n",sep="")
    
  return(invisible(x))
}
