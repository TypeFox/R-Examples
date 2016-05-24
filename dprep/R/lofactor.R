lofactor <-
function (data, k) 
{    if (sum(is.na(data))> 0) 
    stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    data = as.matrix(data)
    distdata = dist.to.knn(data, k)
    p = dim(distdata)[2]
    lrddata = reachability(distdata, k)
    lof = rep(0, p)
    for (i in 1:p) {
        nneigh = distdata[2, i] - distdata[1, i] + 1
        j = seq(0, (nneigh - 1))
        local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
        lof[i] = local.factor
    }
    lof
}
