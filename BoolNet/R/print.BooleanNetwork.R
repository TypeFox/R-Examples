# Custom print function for class BooleanNetwork
print.BooleanNetwork <- function(x, ...)
{
  cat("Boolean network with",length(x$genes),"genes\n\n")
  cat("Involved genes:\n",paste(x$genes,collapse=" "),"\n\n",sep="")
  cat("Transition functions:\n")
  
  mapply(function(gene,interaction)
    {
      # print original expressions read from the files (if available)
      cat(gene," = ",interaction$expression,"\n",sep="")
    },  
    x$genes,x$interactions)
  
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
    
  return(invisible(x))  
}
