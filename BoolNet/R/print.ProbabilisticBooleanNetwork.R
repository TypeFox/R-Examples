print.BooleanNetworkCollection <- function(x, ...)
{
  print.ProbabilisticBooleanNetwork(x, ...)
}

# Custom print function for class ProbabilisticBooleanNetwork
print.ProbabilisticBooleanNetwork <- function(x, ...)
{
  cat("Probabilistic Boolean network with",length(x$genes),"genes\n\n")
  cat("Involved genes:\n",paste(x$genes,collapse=" "),"\n\n",sep="")
  cat("Transition functions:\n")
    
  mapply(function(gene,interaction)
    {
      cat("\nAlternative transition functions for gene ",gene,":\n",sep="")
      # print original expressions read from the files (if available)
      lapply(interaction,function(func)
      {
        cat(gene," = ",func$expression,sep="")
        
        if (!is.null(func$probability) || !is.null(func$error))
        {
          cat(" (")
          
          if (!is.null(func$probability))
          {
            cat(" probability: ",func$probability,sep="")
            if (!is.null(func$error))
              cat(", ")
          }  
          if (!is.null(func$error))
            cat("error: ",func$error,sep="")
          
          cat(")")
        }
        cat("\n")
      })
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
