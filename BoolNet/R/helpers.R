
# Encode a vector of binary values <bin> with <len> bits
# to a decimal number
bin2dec <- function(bin,len)
{
  if (len %% 32 == 0)
    numElts <- as.integer(len / 32)
  else
    numElts <- as.integer(len / 32) + 1

  dec = rep(0,numElts)
  
  dec = .C("bin2dec",as.integer(dec),as.integer(bin),as.integer(len))[[1]]
}

# Decode the <len> low-order bits of <dec> to a vector of binary values,
# where the first entry is the least significant bit
dec2bin <- function(dec,len)
{
  bin = rep(0,len)
  
  bin = .C("dec2bin",as.integer(bin),as.integer(dec),as.integer(len),NAOK=TRUE)[[1]]
}

# Generate a list of all assignments of n variables with N possible values
allcombn <- function(N,n)
{
  rownum = N^n
  sapply(n:1,function(i)
    {
            rep(seq_len(N),each=N^(i-1),len=rownum)
          })
}

# Get DNF representation of a truth table <truthTable>
# using the gene names in <genes>. 
# If <mode> is "canonical", build a canonical DNF.
# If <mode> is "short", join terms to reduce the DNF
getDNF <- function(truthTable, genes, mode = c("short","canonical"))
{
  if (mode[1] == TRUE)
    mode <- (if (length(genes) <= 12) "short" else "canonical")
    
  mode <- match.arg(mode, c("short","canonical"))
  # check for constant functions
  if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
    return("0")
  else
  if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
    return("1")
  
  # generate truth table
  entries <- allcombn(2,length(genes)) - 1
  colnames(entries) <- genes

 
  if (mode == "short")
  {
    # heuristic minimization
    
    # the 1 terms that need to be covered
    uncoveredEntries <- which(truthTable == 1)

    # current conjunction list
    conjunctions <- list()  
    
    while (length(uncoveredEntries) > 0)
    {
      # take an uncovered entry and expand it
      currentEntry <- entries[uncoveredEntries[1],]
          
      for (gene in genes)
      # test for each gene whether it can be eliminated from the term
      {
        geneIdx <- which(names(currentEntry) == gene)
        candidate <- currentEntry[-geneIdx]
        condition <- rep(TRUE,length(truthTable))
        for (i in seq_along(candidate))
        {
          condition <- condition & (entries[,names(candidate)[i]] == candidate[i])
        }
        
        if (length(unique(truthTable[condition])) == 1)
          # eliminate gene
          currentEntry <- currentEntry[-geneIdx]
      }
      
      # determine which truth table result entries are now covered
      eliminatedEntries <- rep(TRUE,length(truthTable))
      for (i in seq_along(currentEntry))
      {
        eliminatedEntries <- eliminatedEntries & 
                             (entries[,names(currentEntry)[i]] == currentEntry[i])
      }
      uncoveredEntries <- setdiff(uncoveredEntries, which(eliminatedEntries))
      
      # remember conjunction
      conjunctions <- c(conjunctions, list(currentEntry))
    }
    return(paste(paste("(",sapply(conjunctions, function(conj)
    {
      paste(mapply(function(gene, val)
      {
        if (val == 1)
          return(gene)
        else
          return(paste("!",gene,sep=""))
      }, names(conj), conj), collapse=" & ")
    }), ")", sep=""), collapse=" | "))
  }
  else
  {
    # canonical DNF
    conjunctions <- apply(entries[truthTable==1,,drop=FALSE],1,function(conj)
    {
      paste("(",paste(sapply(seq_along(conj),function(lit)
      {
        if (conj[lit])
          genes[lit]
        else
          paste("!",genes[lit],sep="")
      }),collapse=" & "),")",sep="")
    })
    return(paste(conjunctions[conjunctions != ""],collapse = " | "))
  }
}

# Retrieves a string representation of an interaction function by either
# building a DNF (if <readableFunction> is false)
# or returning an unspecific function description 
# (if <readableFunction> is true or "canonical" or "short").
# <truthTable> contains the result column of the interaction's truth table, and
# <genes> contains the names of the involved genes.
getInteractionString <- function(readableFunctions,truthTable,genes)
{
  if (readableFunctions != FALSE)
    getDNF(truthTable,genes, readableFunctions)
  else
  {
    if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
      return("0")
    else
    if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
      return("1")
    else
    {
      truthTable <- sapply(truthTable,function(x)
      {
        if (x == 0)
          "0"
        else
        if (x == 1)
          "1"
        else  
          "*"
      })
      paste("<f(",
        paste(genes,collapse=","),"){",
        paste(truthTable,collapse=""),"}>", sep="")
    }
  }    
}

# Reorders a matrix of states <stateMatrix> (with each column being one state)
# in a canonical way such that the "smallest" state is the first.
# This makes attractor representations unique.
canonicalStateOrder <- function(stateMatrix)
{
  smallestIndex <- -1
  smallestVal <- rep(Inf,nrow(stateMatrix))
  for (i in seq_len(ncol(stateMatrix)))
  # iterate over states
  {
    for (j in seq_len(nrow(stateMatrix)))
    # iterate over elements of encoded state
    {
      if (stateMatrix[j,i] < smallestVal[j])
      # determine new minimum
      {
        smallestVal <- stateMatrix[,i]
        smallestIndex <- i
        break
      }
      else
      {
        if (stateMatrix[j,i] > smallestVal[j])
          break
      }
    }
  }
  if (smallestIndex != 1)
    # rearrange matrix
    return(cbind(stateMatrix[,smallestIndex:ncol(stateMatrix),drop=FALSE],
           stateMatrix[,seq_len(smallestIndex-1),drop=FALSE]))
  else
    return(stateMatrix)
}

# Find a child node named <name>
# of <node>, or return NULL if such a node
# does not exist.
# If <throwError>, throw an error if the node does not exist.
xmlFindNode <- function(node,name,throwError=FALSE)
{
    indices <- which(xmlSApply(node, xmlName) ==  name)
    if (length(indices) == 0)
    {
      if (throwError)
        stop(paste("Node \"",name,"\" is required, but missing!", sep=""))
      else
        return(NULL)
    }
    return(xmlChildren(node)[indices])
}

# Remove leading and trailing whitespace characters from <string>
trim <- function(string)
{
  string <- gsub("^[ \t]+", "", string)
  string <- gsub("[ \t]+$", "", string)
  return(string)
}

# Generate an interaction by parsing <expressionString>
# and building the corresponding truth table.
# Here, <genes> is a list of all genes in the network.
# Returns an interaction as used in the BooleanNetwork class.
generateInteraction <- function(expressionString, genes)
{
  res <- .Call("getTruthTable_R", parseBooleanFunction(expressionString, genes), length(genes))
  names(res) <- c("input","func")
  res$expression <- expressionString
  return(res)
}

# Create a DNF from a count predicate 
# ("maj","sumis","sumgt","sumlt").
# <tree> is the expression tree of the predicate
# Returns an expression tree in DNF
expandCountPredicate <- function(tree)
{
  if (tree$operator == "maj")
  {
    if (length(tree$operands) %%2 == 0)
      count <- length(tree$operands) / 2
    else
      count <- floor(length(tree$operands) / 2)
       
    operands <- tree$operands      
  }
  else
  {
    count <- tree$operands[[length(tree$operands)]]$value
    operands <- tree$operands[-length(tree$operands)]
  }
  table <- as.matrix(allcombn(2,length(operands)) - 1)
  tableRes <- as.integer(apply(table,1,function(x)
  {
    switch(tree$operator,
      "maj" = sum(x)  > count,
      "sumgt" = sum(x)  > count,
      "sumlt" = sum(x)  < count,
      "sumis" = sum(x) == count
      )
  }))
  return(parseBooleanFunction(getDNF(tableRes, sapply(operands, stringFromParseTree))))
}

# Get the inputs for an expression tree <tree>.
# If <index> is TRUE, the indices are returned instead of the gene names.
# Returns a vector of gene names.
getInputs <- function(tree, index=FALSE)
{
 res <- switch(tree$type,
    operator = 
    {
      unique(unlist(lapply(tree$operands, getInputs, index=index)))
    },    
    atom = 
    {
      if (index)
        tree$index
      else
        tree$name
    },   
    constant =  
    {
      NULL
    }
  )                     
  return(res) 
}

# Determine the maximum time delays for the genes <genes> 
# in a symbolic expression tree <tree>.
# Returns a vector of time delays for the genes.
maxTimeDelay <- function(tree, genes)
{
  res <- switch(tree$type,
    operator = 
    {
      timeDelays <- sapply(tree$operands, maxTimeDelay, genes=genes)
      
      if (!is.null(dim(timeDelays)))
        apply(timeDelays,1,max)
      else
        timeDelays
    },    
    atom = 
    {
      res <- rep(1,length(genes))
      names(res) <- genes
      res[tree$index] <- -tree$time
      res
    },   
    constant =  
    {
      res <- rep(1,length(genes))
      names(res) <- genes
      res
    }
  )                     
  return(res)  
}

# Convert <geneNames> into valid identifiers
# by replacing special characters
# Returns a vector of adjusted gene names.
adjustGeneNames <- function(geneNames)
{
  # eliminate special characters
  res <- gsub("[^a-zA-Z0-9_]+","_",geneNames)
  # ensure that identifier does not start with a number
  res <- gsub("(^[0-9][a-zA-Z0-9_]*)$","_\\1", res)
  names(res) <- geneNames
  return(res)
}

# Check whether the internal C pointer <ptr> is null
checkNullPointer <- function(ptr)
{
  return(.Call("checkNullPointer",ptr))
}

# Print a synchronous attractor specified by a data frame <attractor>
# The attractor has the index <index> and a basin consisting of <basinSize> states.
# If <activeOnly> is true, states are represented as a list of active genes (otherwise binary vectors).
printSynchronousAttractor <- function(attractor, index, basinSize=NA, activeOnly=FALSE)
{
  if (is.na(basinSize))
    cat("Attractor ",index," is a simple attractor consisting of ", 
        nrow(attractor),
        " state(s)",sep="")
  else
    cat("Attractor ",index," is a simple attractor consisting of ", 
        nrow(attractor),
        " state(s) and has a basin of ", 
        basinSize, " state(s)",sep="")
    
  # print a graphical representation of the attractor cycle
  if (activeOnly)
  {
    cat(".\nActive genes in the attractor state(s):\n")
    
    for (j in seq_len(nrow(attractor)))
    {
      
      state <- paste(colnames(attractor)[which(attractor[j,] == 1)],collapse=", ")
      if (state == "")
        state <- "--"
      cat("State ", j, ": ", state, "\n", sep="")
    }
    cat("\n")
  }
  else
  {
    cat(":\n\n")
    numGenes <- ncol(attractor)
    cat(" |--<", paste(rep("-",numGenes-1), collapse=""),"|\n", sep="")
    cat(" V ", paste(rep(" ",numGenes-1), collapse=""),"  |\n", sep="")
    for (j in seq_len(nrow(attractor)))
    {
      cat(" ", as.integer(attractor[j,]),"   |\n",sep="")
      #cat(" | ", paste(rep(" ",numGenes-1), collapse=""),"  |\n", sep="")
    }
    cat(" V ", paste(rep(" ",numGenes-1), collapse=""), "  |\n", sep="")
    cat(" |-->", paste(rep("-",numGenes-1), collapse=""), "|\n\n", sep="")
  }
  if (!activeOnly)
   cat("\nGenes are encoded in the following order: ",
       paste(colnames(attractor), collapse=" "), "\n\n", sep="")
}

