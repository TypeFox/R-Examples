validation.skewness.kurtosis <-
function(n.NN, skewness.vec=NULL, kurtosis.vec=NULL){
    
    if (missing(n.NN) == TRUE && !is.null(skewness.vec) && !is.null(kurtosis.vec)) {
       stop("Number of continuous variables is not specified !")
    } else
    if ((n.NN < 0) | (floor(n.NN) != n.NN))    {
        stop("Number of continuous variables must be a non-negative integer !")
    } else 
    if (n.NN > 0 && is.null(skewness.vec) && is.null(kurtosis.vec) )         {
        stop("Skewness and kurtosis vectors are not specified while n.NN > 0 !")
    } else
    if (n.NN > 0 && is.null(skewness.vec) && !is.null(kurtosis.vec) )        {
        stop("Skewness vector is not specified while kurtosis.vector is specified and n.NN > 0 !")
    } else
    if (n.NN > 0 && !is.null(skewness.vec) && is.null(kurtosis.vec) )        {
        stop("Kurtosis vector is not specified while skewness.vector is specified and n.NN > 0 !")
    } else
    if (!is.null(skewness.vec) && !is.null(kurtosis.vec)) {
        
         if(n.NN == 0) {
         stop("Skewness and kurtosis vectors are specified while n.NN=0 !")
         } else
         if (n.NN > 0) {
            if( length(skewness.vec)!=length(kurtosis.vec))  {
            stop("Lengths of skewness and kurtosis vectors differ!")
            } else
            if (length(skewness.vec)==length(kurtosis.vec)){
                if( length(skewness.vec)!= n.NN)   {
                 stop("Skewness and kurtosis vectors are misspecified, dimension is wrong!")
                } else
                if (length(skewness.vec)== n.NN)   {
                errorCount= 0         
                for (i in 1:n.NN){
                 if(kurtosis.vec[i] < (skewness.vec[i]^2-2)){
                 cat("\n Kurtosis of continuous variable",i,"must be greater than or equal to",(skewness.vec[i]^2-2),"given its skewness!","\n")
                errorCount = errorCount + 1
                } #for
                } #if
                } #if
                if (errorCount > 0) {
                stop("Range violation occurred in the kurtosis vector!")
                }#if
            }#if
         } #if
     } #if

return(TRUE)

}
