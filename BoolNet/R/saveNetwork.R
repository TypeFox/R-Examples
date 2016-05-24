
# Save <network> to <file> in the BoolNet file format.
# If <generateDNFs> is true or a mode of getDNF(), 
# a new symbolic representation for the interactions
# is generated on the basis of the truth tables (in disjunctive normal form).
# If <generateDNFs> is false, the $expression elements of the interactions are exported.
# If <saveFixed> is true, constant transition functions are exported for fixed genes
# instead of their true transition functions
saveNetwork <- function(network, file, generateDNFs = FALSE, saveFixed = TRUE)
{
  stopifnot(inherits(network,"BooleanNetwork") || inherits(network,"SymbolicBooleanNetwork") ||
            inherits(network,"ProbabilisticBooleanNetwork"))
  expressions <- NULL
  
  if (inherits(network,"BooleanNetwork"))
  {
    if (generateDNFs == FALSE)
    # Check whether all interactions have suitable string representations
    {
      tryCatch(
      {
         # parse the string representations
         lapply(network$interactions,
                function(int)parseBooleanFunction(int$expression))
      },
      error=function(e)
      {
        warning(paste("The transition functions of this network did not contain valid symbolic expressions!",
                "Generating DNF representations from the truth tables!"))
        # There was an error parsing a string representation => generate DNFs              
        generateDNFs <<- TRUE
      })
    }
    
    if (generateDNFs != FALSE)
    {
      expressions <- mapply(function(interaction, gene, fixed)
      {
        if (!saveFixed || fixed == -1)
        {
          table <- allcombn(2, length(interaction$input)) - 1
          expression <- getDNF(interaction$func, 
                               network$genes[interaction$input],
                               generateDNFs)
          return(paste(gene, ", ", expression, sep=""))
        }
        else
          return(paste(gene, ", ", fixed, sep=""))
      }, network$interactions, network$genes, network$fixed)
    }
    else
    {
      expressions <- mapply(function(interaction, gene, fixed)
      {
        if (!saveFixed || fixed == -1)
          return(paste(gene, ", ", interaction$expression, sep=""))
        else
          return(paste(gene, ", ", fixed, sep=""))
      }, network$interactions, network$genes, network$fixed)
    }
    sink(file)
    cat("targets, factors\n")
  }
  else
  if (inherits(network,"SymbolicBooleanNetwork"))
  {
    expressions <- mapply(function(interaction, gene, fixed)
    {
      if (!saveFixed || fixed == -1)
        return(paste(gene, ", ", stringFromParseTree(interaction), sep=""))
      else
        return(paste(gene, ", ", fixed, sep=""))
    }, network$interactions, network$genes, network$fixed)
    sink(file)
    cat("targets, factors\n")
  }
  else
  {
    if (generateDNFs == FALSE)
    # Check whether all interactions have suitable string representations
    {
      tryCatch(
      {
         # parse the string representations
         lapply(network$interactions,
                function(int)lapply(int, function(func)parseBooleanFunction(func$expression)))
      },
      error=function(e)
      {
        warning(paste("The transition functions of this network did not contain valid symbolic expressions!",
                "Generating DNF representations from the truth tables!"))
        # There was an error parsing a string representation => generate DNFs              
        generateDNFs <<- TRUE
      })
    }
    
    if (generateDNFs != FALSE)
    {
      expressions <- mapply(function(interaction, gene, fixed)
      {
        if (!saveFixed || fixed == -1)
        {
          lapply(interaction, function(func)
          {
              table <- allcombn(2, length(func$input)) - 1
              expression <- getDNF(func$func, 
                                   network$genes[func$input],
                                   generateDNFs)
              return(paste(gene, ", ", expression, ", ", func$probability, sep=""))
          })
        }
        else
          return(paste(gene, ", ", fixed, ", 1", sep=""))
      }, network$interactions, network$genes, network$fixed)
    }
    else
    {
      expressions <- mapply(function(interaction, gene, fixed)
      {
        if (!saveFixed || fixed == -1)
        {
          lapply(interaction, function(func)
           {
             return(paste(gene, ", ", func$expression, ", ", func$probability, sep=""))
           })
         }
         else
           return(paste(gene, ", ", fixed, ", 1", sep=""))
         
      }, network$interactions, network$genes, network$fixed)
    }
    expressions <- unlist(expressions)
    sink(file)
    cat("targets, factors, probabilities\n")
  }
  cat(paste(expressions, collapse="\n"),"\n",sep="")
  sink()
}
