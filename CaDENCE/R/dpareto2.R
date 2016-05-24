dpareto2 <-
function(x, scale = 1, shape = 1)
{
    d <- (shape*scale^shape)/((x + scale)^(shape + 1))
    d[x <= 0] <- 0
    d
}
