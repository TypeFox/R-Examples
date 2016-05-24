signorm <-
function (data) 
{
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    classes = data[, d[2]]
    zdata = znorm(data)
    d2 = dim(zdata)
    zdata = zdata[, -d2[2]]
    sigdata = (1 - exp(-zdata))/(1 + exp(-zdata))
    sigdata = cbind(sigdata, classes)
    if (c == "data.frame") 
        sigdata = as.data.frame(sigdata)
    colnames(sigdata) = cnames
    return(sigdata)
}
