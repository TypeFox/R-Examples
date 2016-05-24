acv <-
function (penalty, yy, B, quantile, DD, nb, constmat) 
{
    aa <- asyregpen.lsfit(yy, B, quantile, abs(penalty), DD, 
        nb, constmat)
    score = aa$weight * (yy - B %*% aa$a)^2/(1 - aa$diag.hat.ma)^2
    mean(score[which(is.finite(score))], na.rm = TRUE)
}
