f_rosenbrock <-
function(x){
n <- length(x)
sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}
