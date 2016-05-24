cvgrid <-
function (yy, B, quantile, DD, nb, constmat, types) 
{
    las1 = seq(-2, 4, by = 0.5)
    glatterms = which(types != "parametric")
    lambdas = matrix(las1, nrow = length(las1) * length(glatterms), 
        ncol = 1)
    if (length(glatterms) > 1) 
        for (i in 2:length(glatterms)) lambdas = cbind(lambdas, 
            rep(las1, each = i, times = length(glatterms) - i + 
                1))
    score = rep(0, nrow(lambdas))
    lambdas = 10^lambdas
    penalty = rep(0, length(types))
    for (i in 1:nrow(lambdas)) {
        penalty[glatterms] = lambdas[i, ]
        aa <- asyregpen.lsfit(yy, B, quantile, abs(penalty), 
            DD, nb, constmat)
        score[i] = mean(aa$weight * (yy - B %*% aa$a)^2/(1 - 
            aa$diag.hat.ma)^2, na.rm = T)
    }
    penalty[glatterms] = lambdas[which.min(score), ]
    penalty
}
