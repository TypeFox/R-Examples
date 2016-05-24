## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

kron_tmp <- function(mat1,mat2){
    mat <- mat1 %x% mat2
    return(mat)
}
