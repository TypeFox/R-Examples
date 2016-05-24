kernel_cdf <-
function(x, op_K)
{
    y = pnorm(x, mean = 0, sd = 1)

    res <- y
    return(res)
}
