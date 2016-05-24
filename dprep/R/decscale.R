decscale <-
function (data) 
{
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    classes = data[, d[2]]
    data = data[, -d[2]]
    maxvect = apply(abs(data), 2, max)
    kvector = ceiling(log10(maxvect))
    scalefactor = 10^kvector
    decdata = scale(data, center = FALSE, scale = scalefactor)
    attributes(decdata) = NULL
    decdata = matrix(decdata, dim(data)[1], dim(data)[2])
    decdata = cbind(decdata, classes)
    if (c == "data.frame") 
        decdata = as.data.frame(decdata)
    colnames(decdata) = cnames
    return(decdata)
}
