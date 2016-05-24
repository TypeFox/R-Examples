softmaxnorm <-
function (data) 
{
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    classes = data[, d[2]]
    zdata = znorm(data)
    d2 = dim(zdata)
    zdata = zdata[, -d2[2]]
    softdata = 1/(1 + exp(-zdata))
    softdata = cbind(softdata, classes)
    if (c == "data.frame") 
        softdata = as.data.frame(softdata)
    colnames(softdata) = cnames
    return(softdata)
}
