znorm <-
function (data) 
{
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    classes = data[, d[2]]
    data = data[, -d[2]]
    zdata = scale(data)
    attributes(zdata) = NULL
    zdata = matrix(zdata, dim(data)[1], dim(data)[2])
    zdata = cbind(zdata, classes)
    if (c == "data.frame") 
        zdata = as.data.frame(zdata)
    colnames(zdata) = cnames
    return(zdata)
}
