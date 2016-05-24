## generate uniform distribution random samples
UniformSample <- function(para)
{
    runif(para[1], min=para[2], max=para[3])
}
