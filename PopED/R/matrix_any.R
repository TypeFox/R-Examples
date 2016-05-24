## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

matrix_any<- function(mat){
    return(any(as.logical(mat)))
}
