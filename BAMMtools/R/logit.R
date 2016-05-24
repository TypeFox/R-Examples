logit <- function (x, min = 0, max = 1) 
{
    p <- (x - min)/(max - min)
    log(p/(1 - p))
}
