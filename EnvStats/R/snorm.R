snorm <-
function (q, mean = 0, sd = 1) 
{
    1 - pnorm(q, mean = mean, sd = sd)
}
