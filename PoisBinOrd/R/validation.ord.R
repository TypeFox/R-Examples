validation.ord <-
function(n.O, prop.list = NULL){
  
    if (missing(n.O) == TRUE && !is.null(prop.list)) {
       stop("Number of ordinal variables is not specified!")
    }
    if ((n.O < 0) | (floor(n.O) != n.O))    {
        stop("Number of ordinal variables must be a non-negative integer!")
    } else 
    if (n.O > 0 && is.null(prop.list))        {
        stop("Proportion list is not specified while n.O > 0!")
    } else
    if (!is.null(prop.list)) {
       if(n.O == 0) {
         stop("Proportion list is specified while n.O=0!")
       } else 
       if (n.O > 0 && (length(prop.list) != n.O)) {
        stop("Proportion list is misspecified, dimension is wrong!")
       } else 
       if (n.O > 0 && (length(prop.list) = n.O)) {
        errorCount1= 0
        for (i in 1:n.O){
        for (j in 1:length(prop.list[[i]])){
          if(prop.list[[i]][j] <= 0 | prop.list[[i]][j] >= 1) {
          cat("\n Cumulative proportion for ordinary variable",i,"must be between '0' and '1'!","\n")
          errorCount1 = errorCount1 + 1
          cat("\n")
          } #if
        } #for 
        } #for      
          if (errorCount1 > 0) {
          stop("Range violation occurred in the proportion list!")
          }#if
   
        errorCount2=0
        for (i in 1:n.O){
        if(length(prop.list[[i]])>1)
        for (j in 2:length(prop.list[[i]])){
          if(prop.list[[i]][j] < prop.list[[i]][j-1]) {
          cat("\n Cumulative proportion for ordinary variable",i,"is not in the form of cumulative probabilities!","\n")
          errorCount2 = errorCount2 + 1
          cat("\n")
          } #if
        } #for 
        } #for      
          if (errorCount2 > 0) {
          stop("Range violation occurred in the proportion list!")
          }#if
       } #if

    } #if
return(TRUE)
}
