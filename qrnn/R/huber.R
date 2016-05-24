huber <-
function(x, eps)
{
    h <- ifelse(abs(x)>eps, abs(x)-eps/2, (x^2)/(2*eps))
    h[is.nan(h)] <- 0
    h
}

