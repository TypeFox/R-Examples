# Custom print function for class AttractorInfo
print.SymbolicSimulation <- function(x, activeOnly = FALSE, sequences=FALSE, 
                                     graph=FALSE, attractors=TRUE, ...)
{
  cat("Simulation of a symbolic Boolean network\n")
  
  if (!is.null(x$sequences))
  {
    if (sequences)
    {
      cat(sprintf("Sequences for %d start states:\n", length(x$sequences)))
      
      if (!activeOnly)
        print(x$sequences)
      else
        for (i in seq_along(x$sequences))
        {
          cat("[[",i,"]]\n",sep="")
          seq <- x$sequences[[i]]
          cat(paste(apply(seq,1,function(r)
          {
            if (any(r == 1))
              paste(colnames(seq)[which(r == 1)],collapse=", ")
            else
              "--"
          }), collapse="\n"))
          cat("\n\n")
        }
    }
    else
      cat(sprintf("Sequences for %d start states (print with sequences=TRUE to show them)\n", 
                  length(x$sequences)))
    cat("\n")
  }
  if (!is.null(x$graph))
  {
      if (graph)
    {
      cat(sprintf("Graph containing %d state transitions:\n", nrow(x$graph)))
      print(x$graph, activeOnly=activeOnly, ...)
    }
    else
      cat(sprintf("Graph containing %d state transitions (print with graph=TRUE to show them)\n", 
                  nrow(x$graph)))
    cat("\n")                  
  }
  
  if (!is.null(x$attractors))
  {
    if (attractors)
    {
      cat(sprintf("%d Attractors:\n", length(x$attractors)))
      
      for (i in seq_along(x$attractors))
      {
        if (is.null(x$graph))
          printSynchronousAttractor(x$attractors[[i]], i, activeOnly=activeOnly)
        else
          printSynchronousAttractor(x$attractors[[i]], i, 
                                    sum(x$graph$attractorAssignment == i), activeOnly=activeOnly)
      }
    }      
    else
        cat(sprintf("%d Attractors (print with attractors=TRUE to show them)\n", 
                    length(x$attractors)))  
  }                    
  return(invisible(x))
}
