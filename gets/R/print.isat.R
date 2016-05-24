print.isat <-
function(x, ...)
{
  ##determine type:
  spec <- "mean"

  ##header:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Method: Ordinary Least Squares (OLS)\n")
  if(spec=="mean"){
    cat("Variance-Covariance:", switch(x$aux$vcov.type,
      ordinary = "Ordinary", white = "White (1980)",
      "newey-west" = "Newey and West (1987)"), "\n")
    cat("No. of observations (mean eq.):", x$aux$y.n, "\n")
    cat("Sample (mean eq.):",
      as.character(x$aux$y.index[1]), "to",
      as.character(x$aux$y.index[x$aux$y.n]), "\n")
  }

  ##gum:
  if(spec=="mean"){
    cat("\n")
    cat("GUM mean equation:\n")
    cat("\n")
    print(x$gum.mean)
  }
  if(!is.null(x$gum.variance)){
    cat("\n")
    cat("GUM log-variance equation:\n")
    cat("\n")
    print(x$gum.variance)
  }
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  print(x$gum.diagnostics)

  ##paths:
  cat("\n")
  cat("Paths searched: \n")
  cat("\n")
  if(is.null(x$paths)){
    print(NULL)
  }else{
    for(i in 1:length(x$paths)){
      cat("path",i,":",x$paths[[i]],"\n")
    }
  } #end if(is.null(x$paths))

  ##terminal models and results:
  cat("\n")
  cat("Terminal models: \n")
  if(!is.null(x$terminals)){
    cat("\n")
    for(i in 1:length(x$terminals)){
      cat("spec",i,":",x$terminals[[i]],"\n")
    }
  }
  cat("\n")
  print(x$terminals.results)

  ##specific model:
  if(spec=="mean" && !is.null(x$specific.spec)){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if(!is.null(x$mean.results)){
      print(x$mean.results)
    }
    if(x$specific.spec[1]==0){
      cat("empty\n")
    }
##in the future: use estimate.specific=FALSE more directly?
    if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
      cat("Not estimated\n")
    }
  }
  if(!is.null(x$specific.diagnostics)){
    cat("\n")
    cat("Diagnostics:\n")
    cat("\n")
    print(x$specific.diagnostics)
  }

  ##delete this one? or at least change?:
  ##notes:
#  if(!is.null(x$notes)){
#    cat("\n")
#    cat("Notes:\n")
#    cat("\n")
#    for(i in 1:length(x$notes)){
#      cat("-",x$notes[[i]],"\n")
#    }
#    cat("\n")
#  }
}
