## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

matrix_all<- function(mat){
    return(all(as.logical(mat)))
}
