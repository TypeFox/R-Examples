qpareto2 <-
function(p, scale = 1, shape = 1)
{
    if (any(p <= 0) || any(p >= 1)) stop("bad input for argument 'p'")
    q <- scale*(-1 + (1-p)^(-1/shape))
    q
}
