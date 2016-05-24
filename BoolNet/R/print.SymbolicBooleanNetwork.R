print.SymbolicBooleanNetwork <- function(x,...)
{
  cat("Symbolic representation of a Boolean network\n\n");
  cat("Transition functions:\n")
  for (gene in x$genes)
  {
    cat(gene, " = ", stringFromParseTree(x$interactions[[gene]]), "\n", sep="")
  } 
  
  if (sum(x$fixed != -1) > 0)
  {
    cat("\nKnocked-out and over-expressed genes:\n")
    mapply(function(gene,fixed)
      {
        if (fixed != -1)
          cat(gene," = ",fixed,"\n",sep="")
      },
      x$genes,x$fixed)
  }
    
}

