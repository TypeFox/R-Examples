"orderFeatures" <-
function(x, cl, subset=NULL, FUN=aovFbyrow, values=FALSE){
    if(dim(x)[2]!=length(cl))stop(paste("Dimension 2 of x is",
                  dim(x)[2], "differs from the length of cl (=",
                  length(cl)))
    ## Ensure that cl is a factor & has no redundant levels
    if(is.null(subset))
      cl <- factor(cl)
    else
      cl <- factor(cl[subset])
     if(is.null(subset))
      stat <- FUN(x, cl)
    else
      stat <- FUN(x[, subset], cl)
    ord <- order(-abs(stat))
    if(!values)ord else(list(ord=ord, stat=stat[ord]))
  }

