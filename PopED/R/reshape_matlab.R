## Function written to match MATLAB function
## Author: Andrew Hooker

reshape_matlab <- function(mat,dim1,dim2){
    dim(mat) = c(dim1,dim2)
    return(mat)
}
