`is.toeplitz` <-
function(x){
is.matrix(x)&&(nrow(x)==ncol(x))&&(max(abs(x-toeplitz(x[1,])))<.Machine$double.eps)
}

