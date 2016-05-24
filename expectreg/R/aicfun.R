aicfun <-
function (penalty, yy, B, quantile, DD, nb, constmat) 
{
    aa <- asyregpen.lsfit(yy, B, quantile, abs(penalty), DD, 
        nb, constmat)
    score = length(yy) * log(mean(aa$weight * (yy - B %*% aa$a)^2)) + 
        2 * (1 + sum(aa$diag.hat.ma))
    score
}
