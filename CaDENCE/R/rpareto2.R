rpareto2 <-
function(n, scale = 1, shape = 1)
{
    if (max(length(scale), length(shape)) > 1) 
        stop("parameters must be of length 1")
    scale*(-1 + runif(n)^(-1/shape))
}
