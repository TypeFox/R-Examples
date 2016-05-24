
print.iModel <- function(x, ...){

  cat(sprintf("Model: A %s with %i variables\n", class(x)[1], length(x$varNames)))
  #str(x$varNames)
  ## Model properties
  cat(sprintf(" graphical : %5s  decomposable : %5s\n", x$isGraphical, x$isDecomposable))

  if (x$isFitted){
    dimension <- x$fitinfo$dimension
    #cat("Fit info: \n")
    cat(sprintf(" -2logL    : %14.2f mdim : %4d aic : %12.2f \n",
                -2*x$fitinfo$logL,       dimension["mod.dim"], x$fitinfo$aic))
    cat(sprintf(" ideviance : %14.2f idf  : %4d bic : %12.2f \n",
                x$fitinfo$ideviance,  dimension["idf"], x$fitinfo$bic))
    cat(sprintf(" deviance  : %14.2f df   : %4d \n",
                x$fitinfo$dev,        dimension["df"]))
  }
    
  return(invisible(x))
}

print.dModel <- function(x, ...){

  print.iModel(x)
  
  ## If the model is fitted
  ##
  if (x$isFitted){    
    ## Print warnings about sparsity and adjusments of df's
    ##
    if ( (!x$fitinfo$sparseinfo["df.ok"]) | (!x$fitinfo$sparseinfo["chi2.ok"])) {
      cat("Notice: Table is sparse\n")
      if (!x$fitinfo$sparseinfo["chi2.ok"])
        cat(sprintf("  Asymptotic chi2 distribution may be questionable.\n"))
      
      if (!x$fitinfo$sparseinfo["df.ok"])
        cat(sprintf("  Degrees of freedom can not be trusted.\n"))
      
      if (x$fitinfo$sparseinfo["sparse.df.ok"] & !x$fitinfo$sparseinfo["df.ok"]){
        cat(sprintf("  Model dimension adjusted for sparsity : %d\n",
                    x$fitinfo$dimension["mod.dim.adj"]))
      }
    }
  }
  return(invisible(x))
}













































## ..print.mModel <- function(x,...)
##   {
    
##     cat("Mixed interaction model: \n")
    
##     cat("Model:\n")
##     utils::str(x$glist, give.head=FALSE,no.list=TRUE,comp.str=" ")
    
##     if (x$isFitted){      
##       cat(sprintf("Dimension: %3i df: %3i logL %f -2logL=%f\n",
##                   x$dimension[1], x$dimension[4],x$fitinfo$logL, -2*x$fitinfo$logL))
##     } else {
##       cat(sprintf("Dimension: %3i df: %3i\n", x$dimension[1], x$dimension[4]))
##     }
    
##     ##cat("Object has slots:\n")
##     ##print(names(x))
##     return(invisible(x))
##   }


## print.dModel <- function(x, ...){

##   print.iModel(x)
  
##   ## If the model is fitted
##   ##
##   if (x$isFitted){    
##     ## Print warnings about sparsity and adjusments of df's
##     ##
##     if ( (!x$fitinfo$df.ok) | (!x$fitinfo$chi2.ok)) {
##       cat("Notice: Table is sparse\n")
##       if (!x$fitinfo$chi2.ok)
##         cat(sprintf("  Asymptotic chi2 distribution may be questionable.\n"))
      
##       if (!x$fitinfo$df.ok)
##         cat(sprintf("  Degrees of freedom can not be trusted.\n"))
      
##       if (x$fitinfo$sparse.df.ok & !x$fitinfo$df.ok){
##         cat(sprintf("  Model dimension adjusted for sparsity; mdim : %d\n", x$fitinfo$dim.adj))
##       }
##     }
##   }
##   return(invisible(x))
## }

## print.iModel <- function(x, ...){

##   cat(sprintf("Model: A %s with %i variables\n", class(x)[1], length(x$varNames)))

##   if (x$isFitted){
##     cat("Fit info: \n")
##     cat(sprintf(" -2logL    : %18.8f mdim  : %4d \n",  -2*x$fitinfo$logL,       x$fitinfo$dim.unadj))
##     cat(sprintf(" ideviance : %18.8f idf   : %4d \n",     x$fitinfo$ideviance,  x$fitinfo$idf))
##     cat(sprintf(" deviance  : %18.8f df    : %4d \n",     x$fitinfo$lrt,        x$fitinfo$df))
##     cat(sprintf(" aic       : %14.4f \n",     x$fitinfo$aic))
##     cat(sprintf(" bic       : %14.4f \n",     x$fitinfo$bic))
##   }

##   ## Model properties
##   cat(sprintf("is graphical=%s is decomposable=%s\n", x$isGraphical, x$isDecomposable))

##   return(invisible(x))
## }




  ##   if (x$isFitted){
##     dimension <- x$fitinfo$dimension
##     #cat("Fit info: \n")
##     cat(sprintf(" -2logL    : %14.4f mdim  : %4d \n",
##                 -2*x$fitinfo$logL,       dimension["mod.dim"]))
##     cat(sprintf(" ideviance : %14.4f idf   : %4d \n",
##                 x$fitinfo$ideviance,  dimension["idf"]))
##     cat(sprintf(" deviance  : %14.4f df    : %4d \n",
##                 x$fitinfo$lrt,        dimension["df"]))
##     cat(sprintf(" aic       : %14.4f \n bic       : %14.4f \n",
##                 x$fitinfo$aic,        x$fitinfo$bic))
##   }

  
##   if (x$isFitted){
##     dimension <- x$fitinfo$dimension
##     #cat("Fit info: \n")
##     cat(sprintf(" -2logL    : %14.4f ideviance : %14.4f deviance  : %14.4f\n",
##                 -2*x$fitinfo$logL, x$fitinfo$ideviance, x$fitinfo$lrt))
##     cat(sprintf(" mdim      : %4d   idf : %4d  df : %4d\n",
##                 dimension["mod.dim"], dimension["idf"], dimension["df"] ))
    
## ##     cat(sprintf("  df    : %4d \n",     x$fitinfo$lrt,        dimension["df"]))
## ##     cat(sprintf(" aic       : %14.4f \n bic       : %14.4f \n",  x$fitinfo$aic,        x$fitinfo$bic))
##   }
