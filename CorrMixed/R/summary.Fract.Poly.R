summary.Fract.Poly <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  x <- Object
  
  if (is.null(Object$Results.M1)==FALSE){
    cat("Best fitting model for m=1: \n")
    print(Object$Results.M1[order(Object$Results.M1[,2]),][1:1,])
  }
  
  if (is.null(Object$Results.M2)==FALSE){
    cat("\nBest fitting model for m=2: \n")
    print(Object$Results.M2[order(Object$Results.M2[,3]),][1:1,])
  }
  
  if (is.null(Object$Results.M3)==FALSE){
    cat("\nBest fitting model for m=3: \n")
    print(Object$Results.M3[order(Object$Results.M3[,4]),][1:1,])
  }
  
  if (is.null(Object$Results.M4)==FALSE){
    cat("\nBest fitting model for m=4: \n")
    print(Object$Results.M4[order(Object$Results.M4[,5]),][1:1,])
  }
  
  if (is.null(Object$Results.M5)==FALSE){
    cat("\nBest fitting model for m=5: \n")
    print(Object$Results.M5[order(Object$Results.M5[,6]),][1:1,])
  }
}  


