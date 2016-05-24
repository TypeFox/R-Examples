validation.bin <-
function(n.BB, prop.vec = NULL){
  
    if (missing(n.BB) == TRUE && !is.null(prop.vec)) {
       stop("Number of binary variables is not specified !")
    }
    if ((n.BB < 0) | (floor(n.BB) != n.BB))    {
        stop("Number of binary variables must be a non-negative integer !")
    } else 
    if (n.BB > 0 && is.null(prop.vec))         {
        stop("Proportion vector is not specified while n.BB > 0 !")
    } else
    if (!is.null(prop.vec)) {
       if(n.BB == 0) {
         stop("Proportion vector is specified while n.BB=0")
       } else 
       if (n.BB > 0 && (length(prop.vec) != n.BB)) {
        stop("Proportion vector is misspecified, dimension is wrong!")
       } else 
       if (n.BB > 0 && (length(prop.vec) = n.BB)) {
        errorCount= 0
        for (i in 1:n.BB){
          if(prop.vec[i] <= 0 | prop.vec[i] >= 1) {
          cat("\n Proportion for binary variable",i,"must be between '0' and '1!!","\n")
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
