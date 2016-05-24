print.ContaminatedMixt <- function(x,...){
  if (length(x$models) >0) {
    best <- whichBest(x,...)
    best.unique <- unique(best)  
    for (i in seq_len(length(best.unique))){
      if (length(x$models)>1){
        b <- best==best.unique[i]        
        m <- paste(names(best)[b],collapse=", ")
        m <- paste("\nBest model according to", m, "has")
      } else m <- "\nEstimated one model with"
      m <- paste(m, "G =", x$models[[best.unique[i]]]$G,"group(s)")
      if (!is.null(x$models[[best.unique[[i]]]]$model)) 
        m <- paste(m, "and parsimonious structure",  x$models[[best.unique[[i]]]]$model)
           
      cat(m,"\n")
    }
  }
  else cat("No models have been estimated.")
  invisible(x)
}

