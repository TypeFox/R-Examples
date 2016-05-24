##
## INPUT:
## x: numeric vector of values in [0, 1]
##
## RETURN:
## numeric vector, ensuring all values of x are numerically inside (0, 1)
##
finiteProbs <- function(x)
{
    x <- replace(x, x < .Machine$double.eps, .Machine$double.eps)
    x <- replace(x, x > 1 - .Machine$double.neg.eps,
                 1 - .Machine$double.neg.eps)
    return(x)
}
