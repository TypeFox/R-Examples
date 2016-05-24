"loo" =
function (X, Y, lambda1, lambda2) 
{
    n = nrow(X)
    xscore = vector(mode = "numeric", length = n)
    yscore = vector(mode = "numeric", length = n)
    for (i in 1:n) {
        Xcv = X[-i, ]
        Ycv = Y[-i, ]
        res = rcc(Xcv, Ycv, lambda1, lambda2)
        xscore[i] = t(X[i, ]) %*% res$xcoef[, 1]
        yscore[i] = t(Y[i, ]) %*% res$ycoef[, 1]
    }
    cv = cor(xscore, yscore, use = "pairwise")
    return(invisible(cv))
}

