summary.PCAmix <-
  function(object, ...)
  {
    x <- object
    if (!inherits(x, "PCAmix")) 
      stop("use only with \"PCAmix\" objects")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    n <- x$rec$n
    p1 <- x$rec$p1
    p <- x$rec$p
    p2 <- p-p1
    if (colnames(x$sqload)[1]=="dim1.rot")  cat("Method = rotation after ") else  cat("Method = ")
    if (p1==p) cat("Principal Component Analysis (PCA)") 
    if (p1==0) cat("Multiple Correspondence Analysis (MCA)") 
    if ((p1!=0) && (p1!=p)) cat("Factor Analysis of mixed data (FAmix)") 
    cat("\n")
    cat("\n")
    cat("Data:", "\n")
    cat(paste("   number of observations: ",n),sep=" ") 
    cat("\n")
    if 	((p1!=0)&& (p2==0)) {
      cat(paste("   number of variables: ",p1),sep=" ")  
      cat("\n")
    }
    if 	((p1==0)&& (p2!=0)) {
      cat(paste("   number of variables: ",p2),sep=" ")  
      cat("\n")
    }
    if 	((p1!=0)&& (p2!=0)) {
      cat(paste("   number of  variables: ",p),sep=" ")
      cat("\n")
      cat(paste("        number of numerical variables: ",p1),sep=" ")   
      cat("\n")
      cat(paste("        number of categorical variables: ",p2),sep=" ")   
      cat("\n")
    }
    cat("\n")
    if (colnames(x$sqload)[1]=="dim1.rot") cat("Squared loadings after rotation:")
    else cat("Squared loadings :")
    cat("\n")
    print(round(x$sqload,digits=2))
    cat("\n")
    cat("\n")
    
  }

