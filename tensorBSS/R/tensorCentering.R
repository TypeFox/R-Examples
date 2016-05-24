tensorCentering <-
function(x){
  r <- length(dim(x)) - 1
  sweep(x, 1:r, apply(x, 1:r, mean), '-')
}
