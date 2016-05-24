as.matrix.mvabund <- function(x, ...){

    if(!is.mvabund(x))
    stop("The function 'as.matrix.mvabund' can only be used for a mvabund object.")

    if(is.null(dim(x))) x <- c(x)

   	x <- as.matrix(unabund(x))
    return(x)

}
