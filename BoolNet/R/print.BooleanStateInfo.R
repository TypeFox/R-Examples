# Custom print function for class BooleanStateInfo
print.BooleanStateInfo <- function(x, activeOnly=FALSE, ...)
{
  # Create a TransitionTable object and print it
  cat("Transition table of Boolean network\n\n")
   
  transitionTable <- .internal.getTransitionTable(x)
  
  print(transitionTable, activeOnly=activeOnly, ...)
}
