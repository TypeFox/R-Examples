##
## Methods for showing S4-class objects
## ====================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod(f = "show", signature(object = "GoGARCH"), definition = function(object){
  title <- "*** GO-GARCH ***"
  stars <- paste(rep("*", nchar(title)), collapse = "")
  cat("\n")
  cat(paste(stars, "\n"))
  cat(paste(title, "\n"))
  cat(paste(stars, "\n"))  
  cat("\n")
  cat(paste("Components estimated by:", object@estby))
  cat("\n")
  cat(paste("Dimension of data matrix:", paste("(", nrow(object@X), " x ", ncol(object@X), ").", sep = "")))
  cat("\n")
  cat(paste("Formula for component GARCH models:", paste(as.character(object@garchf), collapse = " "), "\n"))
  cat("\n")  
  if(length(object@U) != 0){
    cat("Orthogonal Matrix U:\n")
    print(object@U)
    cat("\n")
    cat("Linear Map Z:\n")
    print(object@Z)
    cat("\n")
  }
  cat("Estimated GARCH coefficients:\n")   
  print(coef(object))
  cat("\n")
  cat("Convergence codes of component GARCH models:\n")
  print(converged(object))
  invisible(object)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "show", signature(object = "Goestica"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "show", signature(object = "Goestmm"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "show", signature(object = "Goestnls"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "show", signature(object = "Goestml"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Orthom"
##
setMethod(f = "show", signature(object = "Orthom"), function(object) print(object@M))
##
## Method definition for objects of class "Gosum"
##
setMethod(f = "show", signature(object = "Gosum"), definition = function(object){
  title <- "*** Summary of GO-GARCH Model ***"
  stars <- paste(rep("*", nchar(title)), collapse = "")
  cat("\n")
  cat(paste(stars, "\n"))
  cat(paste(title, "\n"))
  cat(paste(stars, "\n"))
  cat("\n")
  cat(paste("Used object:", object@name))
  cat("\n")  
  cat("\n")
  cat(paste("Components estimated by:", object@method))
  cat("\n")
  cat(paste("Formula for component GARCH models:", paste(as.character(object@model), collapse = " "), "\n"))
  cat("\n")  
  if(length(object@Zinv) != 0){
    cat("The Inverse of the Linear Map Z:\n")
    print(object@Zinv)
    cat("\n")
  }
  cat("\n")
  title2 <- "*** Estimated Component GARCH models ***"
  stars2 <- paste(rep("*", nchar(title2)), collapse = "")
  cat("\n")
  cat(paste(stars2, "\n"))
  cat(paste(title2, "\n"))
  cat(paste(stars2, "\n"))
  cnames <- names(object@garchc)
  n <- length(object@garchc)
  for(i in 1:n){
    cat("\n")
    cat(paste(cnames[i], "\n"))
    print(object@garchc[[i]])  
  }
  invisible(object)
})
##
## Method definition for objects of class "Gopredict"
##
setMethod(f = "show", signature = "Gopredict", function(object){
  title <- "*** Forecasts of GO-GARCH Model ***"
  stars <- paste(rep("*", nchar(title)), collapse = "")
  cat("\n")
  cat(paste(stars, "\n"))
  cat(paste(title, "\n"))
  cat(paste(stars, "\n"))  
  cat("\n")
  if(nrow(object@Xf) <= 10){
    cat("Conditional variances:\n")
    print(cvar(object))
    cat("\n")
    cat("\n")
    cat("Forecasts of Mean Equation:\n")
    print(object@Xf)
    cat("\n")    
  } else {
    cat("Head of conditional variances:\n")
    print(head(cvar(object)))
    cat("\n")
    cat("\n")
    cat("Forecasts of Mean Equation:\n")
    print(head(object@Xf))
    cat("\n")    
  } 
})
##
## Method definition for objects of class "Goinit"
##
setMethod(f = "show", signature = "Goinit", function(object){
  title <- "*** Object of class Goinit ***"
  stars <- paste(rep("*", nchar(title)), collapse = "")
  cat("\n")
  cat(paste(stars, "\n"))
  cat(paste(title, "\n"))
  cat(paste(stars, "\n"))  
  cat("\n")
  if(length(object@X) != 0){
     cat("Head of data matrix X:\n")
     print(head(object@X), quote = FALSE)
     cat("\n")
  } 
  if(length(object@P) != 0){
    cat("Projection matrix P:\n")
    print(object@P, quote = FALSE)
    cat("\n")
  }
  if(length(object@Dsqr) != 0){
    cat("Square root of eigenvalues Dsqr:\n")
    print(formatC(object@Dsqr), quote = FALSE)
    cat("\n")
  }
  cat(paste("Formula for component GARCH models:\n", paste(as.character(object@garchf), collapse = " "), "\n"))
}
)

