fun.model.lm.optim <-
function (value, x, y, param) 
{
    len <- length(value)
    gld.val <- fun.mean.convert(c(0, value[(len - 2):len]), param)
    resid <- y - x %*% value[-c((len - 2):len)]
    resid <- resid - mean(resid)
    r <- sum(log(dgl(resid, gld.val, param = param))) * -1
    return(r)
}
