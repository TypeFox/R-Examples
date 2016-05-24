# Create LaTeX state tables of all attractors in <attractorInfo>.
# Genes are grouped according to <grouping>.
# An additional title can be supplied in <title>.
# If <plotFixed> is set, fixed variables are included in the plot.
# <onColor> and <offColor> specify the colors of ON/1 and OFF/0 states.
# <file> is the name of the output LaTeX document.
attractorsToLaTeX <- function (attractorInfo, subset, title = "", grouping = list(), plotFixed = TRUE, 
        onColor="[gray]{0.9}",offColor="[gray]{0.6}", 
        reverse=FALSE, file="attractors.tex")  
{
  stopifnot(inherits(attractorInfo,"AttractorInfo") || inherits(attractorInfo, "SymbolicSimulation"))
  
  if (inherits(attractorInfo,"AttractorInfo"))
  {
    numGenes <- length(attractorInfo$stateInfo$genes)
    geneNames <- attractorInfo$stateInfo$genes
  }
  else
  {
    numGenes <- ncol(attractorInfo$attractors[[1]])
    geneNames <- colnames(attractorInfo$attractors[[1]])
  } 
  
  if (missing(subset))
      subset <- seq_along(attractorInfo$attractors)
  else
    if (any(subset > length(attractorInfo$attractors)))
      stop("You specified an attractor index that is greater than the total number of attractors in 'subset'!")

  # escape "_" in LaTeX
  genes = gsub("_", "\\_", attractorInfo$stateInfo$genes)
  
  # determine list of genes to be plotted
  whichFixed <- which(attractorInfo$stateInfo$fixedGenes != -1)
  if (plotFixed | (length(whichFixed) == 0))
    plotIndices <- seq_len(numGenes)
  else
    plotIndices <- (seq_len(numGenes))[-whichFixed]
  
  if (inherits(attractorInfo,"AttractorInfo"))
  {
    # convert decimal state numbers to binary state matrices (one for each attractor)
    binMatrices <- lapply(attractorInfo$attractors,function(attractor)
            {
              res <- matrix(apply(attractor$involvedStates,2,function(state)
                dec2bin(state,numGenes)[plotIndices]),nrow=length(plotIndices))
            })

    # count the numbers of attractors with equal lengths
    attractorLengths <- sapply(attractorInfo$attractors,function(attractor)
                               {
                                  if (is.null(attractor$initialStates))
                                  # simple attractor
                                    ncol(attractor$involvedStates)
                                  else
                                  # complex attractor => extra treatment
                                    -1
                               })
  }
  else
  {
    binMatrices <- lapply(attractorInfo$attractors, t)
    attractorLengths <- sapply(binMatrices, ncol)
  }
  
  lengthTable <- table(attractorLengths)
  lengthTable <- lengthTable[as.integer(names(lengthTable)) != -1]
  
  # Open output file, and print header
  sink(file)
  cat("% Please include packages tabularx and colortbl in your master document:\n",
      "% \\usepackage{tabularx,colortbl}\n\n\n",sep="")
      
  res <- lapply(seq_along(lengthTable),function(i)
  # accumulate all attractors with equal length in one matrix and plot them
  {
     len <- as.integer(names(lengthTable)[i])
     attractorIndices <- intersect(which(attractorLengths == len), subset)
     if (length(attractorIndices) > 0)
     {
      # build accumulated matrix     
      totalMatrix <- c()
      for (mat in binMatrices[attractorIndices])
      {
        totalMatrix <- cbind(totalMatrix,mat)
      }
      rownames(totalMatrix) <- geneNames[plotIndices]
      colnames(totalMatrix) <- sapply(attractorIndices,function(i)paste("Attr",i,".",seq_len(len),sep=""))
    
      if(length(grouping)>0)
      {
         # reorder genes according to the supplied groups
        totalMatrix <- totalMatrix[unlist(grouping$index),]
        separationPositions <- c(1,cumsum(sapply(grouping$index,length)+1))
      }
      else
        separationPositions <- c()
    
      # output table header
      cat("\\begin{table}[ht]\n",
           "\\begin{center}\n",
           "\\caption{",
           title, "Attractors with ",len," state(s)}\n",
           "\\begin{tabularx}{\\linewidth}{l", 
         paste(rep(paste(rep(">{\\centering\\arraybackslash}X", 
              len), collapse = " "),length(intersect(which(attractorLengths == len),subset))),collapse="|"), 
      "}\\hline\n",
           sep="")
    
      cat("\t&\t",paste(paste("\\multicolumn{",len,"}{c}{Attr. ",intersect(which(attractorLengths == len),subset),"}",
            sep=""),collapse="\t&\t"),"\\\\\n")    
    
      # output active and inactive states
      if (reverse)
        indices <- rev(seq_len(nrow(totalMatrix)))
      else
        indices <- seq_len(nrow(totalMatrix))
  
      for(j in indices)
      {
        separator <- which(separationPositions==j)
        if (length(separator) != 0)
        {
          cat("\\hline \\multicolumn{",ncol(totalMatrix) + 1,"}{c}{",grouping$class[separator],"}\\\\ \\hline \n",sep="")
        }
        cat("\\textbf{",rownames(totalMatrix)[j],"}\t&\t",sep="")
        for(i in seq_len(ncol(totalMatrix)))
        {
          if(totalMatrix[j,i] == 1)
            cat("\\cellcolor",onColor,"1",sep="")
          else
            cat("\\cellcolor",offColor,"0",sep="")
          if (i < ncol(totalMatrix))
            cat("\t&\t")
        }
        cat("\\\\\n")
      }
    
      # output frequency of attractor (basin size / number of states)
      if (inherits(attractorInfo,"AttractorInfo"))
      {
        if (is.null(attractorInfo$stateInfo$table))
            freq <- rep(NA, length(attractorIndices))
        else
          freq <- round(sapply(attractorInfo$attractors[attractorIndices],
              function(attractor)attractor$basinSize/ncol(attractorInfo$stateInfo$table)) * 100,2)
      }
      else
      {
        if (!is.null(attractorInfo$graph))
        {
          freq <- round(sapply(attractorIndices, 
                        function(i)sum(attractorInfo$graph$attractorAssignment == i)/
                        nrow(attractorInfo$graph)) * 100,2)
        }
        else
          freq <- rep(NA, length(attractorIndices))
      }

      if (!isTRUE(all(is.na(freq))))
      {
        cat("\\hline Freq.\t&\t",paste(paste("\\multicolumn{",len,"}{c}{",freq,"\\%}",
              sep=""),collapse="\t&\t"),"\\\\\n")
      }

      cat("\\hline\\end{tabularx}\n\\end{center}\n",
          "\\end{table}\n\n",sep="")

      totalMatrix
    }
  })
  
  # return a list of accumulated matrices
  sink()
  names(res) <- names(lengthTable)
  return(res)
}
