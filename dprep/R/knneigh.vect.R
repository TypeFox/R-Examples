knneigh.vect <-
function (x, data, k) 
{
    temp = as.matrix(data)
    numrow = dim(data)[1]
    dimnames(temp) = NULL
    difference <- scale(temp, x, FALSE)
    dtemp <- drop(difference^2 %*% rep(1, ncol(data)))
    dtemp = sqrt(dtemp)
    order.dist <- order(dtemp)
    nndist = dtemp[order.dist]
    knndist = nndist[k + 1]
    neighborhood = drop(nndist[nndist <= knndist])
    neighborhood = neighborhood[-1]
    numneigh = length(neighborhood)
    index.neigh = order.dist[1:numneigh + 1]
    num1 = length(index.neigh) + 3
    num2 = length(index.neigh) + numneigh + 2
    neigh.dist = c(num1, num2, index.neigh, neighborhood)
    return(neigh.dist)
}
