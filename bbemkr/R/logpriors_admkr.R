logpriors_admkr = function(h2, data_x)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    logf = 0
    hn = data_num^(-1/(4 + dm))
    bn = data_num^(-3/(2 * dm + 11))
    dm = length(h2) - 1
    for(i in 1:dm)
    {
        logf = logf + logpriorh2(h2[i]*hn^2) + log(hn^2)
    }
    logf = logf + logpriorh2(h2[dm+1]*bn^2) + log(bn^2)
    return(logf)
}
