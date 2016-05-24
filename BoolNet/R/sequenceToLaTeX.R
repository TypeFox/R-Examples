
sequenceToLaTeX <- function(network, startState, includeAttractorStates = c("all","first","none"),
                            sequence, title="", grouping = list(), plotFixed = TRUE,
                            onColor="[gray]{0.9}",offColor="[gray]{0.6}", highlightAttractor=TRUE,
                            reverse=FALSE, 
                            file="sequence.tex")
{
  if (!missing(network))
  {
    stopifnot(inherits(network,"BooleanNetwork") || inherits(network,"SymbolicBooleanNetwork"))
    if (missing(startState) || !missing(sequence))
      stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
    
    sequence <- getPathToAttractor(network = network,
                                   state = startState, 
                                   includeAttractorStates = includeAttractorStates)

    numGenes <- length(network$genes)                                   
    whichFixed <- which(network$fixed != -1)
    if (plotFixed | (length(whichFixed) == 0))
      plotIndices <- seq_len(numGenes)
    else
      plotIndices <- seq_len(numGenes)[-whichFixed]
    
    attractor <- attributes(sequence)$attractor  
    sequence <- sequence[,plotIndices]
    attributes(sequence)$attractor <- attractor                                  
  }
  else
  {
    if (missing(sequence) || !missing(startState))
        stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
    
    attractor <- attributes(sequence)$attractor          
  }
  
  # escape "_" in LaTeX
  genes = gsub("_", "\\_", colnames(sequence))
  
  # determine list of genes to be plotted
  
  # Open output file, and print header
  sink(file)
  cat("% Please include packages tabularx and colortbl in your master document:\n",
      "% \\usepackage{tabularx,colortbl}\n\n\n",sep="")
      
  totalMatrix <- t(sequence)
  colnames(totalMatrix) <- seq_len(ncol(totalMatrix))
  
  if(length(grouping)>0)
  {
     # reorder genes according to the supplied groups
    totalMatrix <- totalMatrix[unlist(grouping$index),]
    separationPositions <- c(1,cumsum(sapply(grouping$index,length)+1))
  }
  else
    separationPositions <- c()

  if (highlightAttractor && !is.null(attractor))
  {
      header <- paste(paste(rep(">{\\centering\\arraybackslash}X", 
                            min(attractor) - 1), collapse = " "),
                      "|",
                      paste(rep(">{\\centering\\arraybackslash}X", 
                            length(attractor)), collapse = " "),
                      "|", sep="")
  }
  else
    header <- paste(rep(">{\\centering\\arraybackslash}X", 
          ncol(totalMatrix), collapse = " "))
  # output table header
  cat("\\begin{table}[ht]\n",
       "\\begin{center}\n",
       "\\caption{",title,"}\n",
       "\\begin{tabularx}{\\linewidth}{l", 
       header, 
  "}\\hline\n", sep="")
       
  cat("\\textbf{Time}\t&\t",paste(seq_len(ncol(totalMatrix)),collapse="\t&\t"),"\\\\")
  
   if(length(grouping) == 0)
     cat("\\hline\n")
   else
     cat("\n")  

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
      cat("\\hline \\multicolumn{",ncol(totalMatrix) + 1,
          "}{c}{",grouping$class[separator],"}\\\\ \\hline \n",sep="")
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
  
  cat("\\hline")
  if (highlightAttractor && !is.null(attractor))
  {    
    cat("\\multicolumn{",min(attractor),
        "}{c|}{}\t&\t\\multicolumn{",length(attractor),
        "}{c|}{Attractor}\\\\\\cline{",min(attractor) + 1, "-", max(attractor) + 1, "}\n",sep="")
  }

  cat("\\end{tabularx}\n\\end{center}\n",
      "\\end{table}\n\n",sep="")

  sink()
  # return the matrix  
  return(totalMatrix)
}

