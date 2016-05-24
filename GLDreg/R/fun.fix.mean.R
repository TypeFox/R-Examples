fun.fix.mean <-
function (value, data, param) 
{
    gld.val <- fun.mean.convert(c(0, value), param = param, val = mean(data))
    r <- sum(log(dgl(data, gld.val, param = param))) * -1
    return(r)
}
