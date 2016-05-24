hnorm <-
function (x, mean = 0, sd = 1) 
{
    dnorm(x, mean = mean, sd = sd)/snorm(x, mean = mean, sd = sd)
}
