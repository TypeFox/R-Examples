"isocir" <- function(cirmeans =NULL, SCE = NULL, CIRE = NULL, pvalue = NULL, kappa = NULL){
   result <- list()
   result$CIRE <- CIRE
   result$cirmeans <- cirmeans
   result$SCE <- SCE
   if (!is.null(pvalue)){result$pvalue <- pvalue}
   if (!is.null(kappa)){result$kappa <- kappa}
   class(result) <- "isocir"
   return(result)
}

"is.isocir" <- function(x){inherits(x, "isocir")}

"print.isocir" <- function(x, decCIRE = 3, decpvalue = 4, deckappa = 2, ...){
   pvalue <- NULL
   kappa <- NULL
   CIRE <- x$CIRE
   cat("\n Circular Isotonic Regression Estimator (CIRE):\n")
   for (i in 1:length(CIRE)){
    cat("   ", round(CIRE[[i]], decCIRE),"\n")  
   }
   cat("\n Sum of Circular Errors: SCE = ", round(x$SCE, decCIRE), "\n")
   if(!is.null(x$cirmeans)){
     cat("Invisible: unrestricted circular means; these can be\n")
     cat("obtained via $cirmeans\n")
     }
   cat("\n")
   if(!is.null(x$pvalue)){pvalue <- round(x$pvalue, decpvalue)}
   if(!is.null(x$kappa)){kappa <- round(x$kappa, deckappa)}
   if(!all(is.null(pvalue), is.null(kappa))){
    cat("pvalue = ", pvalue)    
    cat("\n kappa = ", kappa, "\n")
    cat(attr(x,"estkappa"), "\n\n")
   }
   #NextMethod("print", x, ...)
   invisible(x)
}