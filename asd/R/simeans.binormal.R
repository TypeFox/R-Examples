simeans.binormal <-
function (n = n, means = means, vars = vars, corr = corr) 
{
    var.mat <- diag(vars)/n
    sd.mat <- sqrt(var.mat)
    cor.mat <- matrix(1, ncol = 2, nrow = 2)
    cor.mat[1, 2] <- corr
    cor.mat[2, 1] <- cor.mat[1, 2]
    varcov.mat <- sd.mat %*% cor.mat %*% sd.mat
    samps <- rmvnorm(n = 1, mean = means, sigma = varcov.mat)
    list(samp1 = samps[, 1], samp2 = samps[, 2])
}
