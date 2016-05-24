#' @export
print.frfast <- function(x = model, ...) {
  
  if (inherits(x, "frfast")) {
    model <- x
    cat("\nCall:\n")
    print(model$call)
    cat("", "\n")
    
    if (model$nmodel == 1) {
      m <- "Nonparametric"
    } else {
      m <- "Allometric"
    }
    cat("*********************************************", "\n")
    cat(m, "Model", "\n")
    cat("*********************************************", "\n")
    cat("\nNumber of Observations: ")
    cat(format(model$n))
    cat("\n")
    
    if (length(unique(model$fmod)) != 1) {
      cat("\nNumber of Factors: ")
      cat(format(length(unique(model$fmod))))
      cat("\n")
    }
    
    # cat('\nResidual Standar Error: ') cat('\n')
    
    cat("\nNumber of Bootstrap Repeats: ")
    cat(format(model$nboot))
    cat("\n")
    
    if (model$nmodel == 1) {
      cat("\nType of Nonparametric Smoother: ")
      cat(format(model$smooth))
      cat("\n")
    }
      
    
    if (model$nmodel == 1 & model$smooth == "kernel") {
      if (model$nf != 1) {
        cat("\nBandwidth: ")
        cat(format(c(model$h0, model$h), digits = 2))
        cat("\n")
      } else {
        cat("\nBandwidth: ")
        cat(format(c(model$h0), digits = 2))
        cat("\n")
      }
      
      cat("\nKernel Function: ")
      if (model$kernel == 1) 
        cat("Epanechnikov")
      if (model$kernel == 2) 
        cat("Triangular")
      if (model$kernel == 3) 
        cat("Gaussian")
    }
    
    
    
    
    
    
    
    
    
    
  }else{
    stop("Argument x must be either frfast object.")
  }
}