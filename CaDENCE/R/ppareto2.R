ppareto2 <-
function(q, scale = 1, shape = 1)
{
    p <- 1 - (1 + q/scale)^(-shape)
    p[q <= 0] <- 0
    p
}
