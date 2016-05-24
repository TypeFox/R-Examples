validation.bin <-
function(n.B, prop.vec = NULL){
  
    if (missing(n.B) == TRUE && !is.null(prop.vec)) {
       stop("Number of binary variables is not specified!")
    }
    if ((n.B < 0) | (floor(n.B) != n.B))         {
        stop("Number of binary variables must be a non-negative integer!")
    } else 
    if (n.B > 0 && is.null(prop.vec))            {
        stop("Proportion vector is not specified while n.B > 0!")
    } else
    if (!is.null(prop.vec)) {
       if(n.B == 0) {
         stop("Proportion vector is specified while n.B=0!")
       } else 
       if (n.B > 0 && (length(prop.vec) != n.B))  {
        stop("Proportion vector is misspecified, dimension is wrong!")
       } else 
       if (n.B > 0 && (length(prop.vec) = n.B))   {
        errorCount= 0
        for (i in 1:n.B){
          if(prop.vec[i] <= 0 | prop.vec[i] >= 1) {
          cat("\n Proportion for binary variable",i,"must be between '0' and '1'!","\n")
          errorCount = errorCount + 1
          cat("\n")
          } #if
        } #for      
          if (errorCount > 0) {
         stop("Range violation occurred in the proportion vector!")
         }#if
       } #if

    } #if
return(TRUE)
}
