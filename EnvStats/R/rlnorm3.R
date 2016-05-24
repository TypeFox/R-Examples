rlnorm3 <-
function (n, meanlog = 0, sdlog = 1, threshold = 0) 
{
    rlnorm(n = n, meanlog = meanlog, sdlog = sdlog) + threshold
}
