vvalen1 <-
function (data, classn) 
{
    p = dim(data)[2]
    data1 = data[data[, p] == classn, ]
    data1 = znorm(data1)
    med1 = apply(data1, 2, median)
    medm = as.matrix(rep(1, p)) %*% med1
    data2 = data1 - medm
    data2 = data2[, 1:(p - 1)]
    means = sqrt(apply(data2^2, 1, sum))
    return(means)
}
