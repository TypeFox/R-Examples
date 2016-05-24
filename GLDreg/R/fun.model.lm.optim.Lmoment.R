fun.model.lm.optim.Lmoment <-
function (value, x, y, param) 
{
    len <- length(value)
    gld.val <- fun.mean.convert(c(0, value[(len - 2):len]), param)
    resid <- y - x %*% value[-c((len - 2):len)]
    resid <- resid - mean(resid)
    r<-sum((Lmoments(resid) - fun.lm.theo.gld(gld.val[1], 
            gld.val[2], gld.val[3], gld.val[4], param))^2)
    return(r)
}
